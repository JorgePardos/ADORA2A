"""End-to-end reproducible execution pipeline for the ADORA2A Retrospective Audit (v2).

Runs:
  - Phase 1: Structural audit across all 188 ADORA2A PDB structures, pocket residue checks,
             and objective crystallographic water audit (`inspect_pocket_waters` on 4EIY)
  - Phase 2: ChEMBL251 curation with provenance, PAINS/Brenk contextual audit, & 7D decoy gate
  - Phase 3: Clean A/B experiment on 4EIY dry vs 4EIY hydrated:
             * 10-seed cognate ZM241385 redocking (symmetry-corrected RMSD & ΔE distribution)
             * 1,000-compound retrospective benchmark across 5 series (1D, 2D ECFP4 + Murcko split,
               3D Shape + USRCAT, AutoDock Vina 4EIY dry, AutoDock Vina 4EIY hydrated)
             * 60-ligand experimental pKi correlation panel across >5 orders of magnitude

Usage:
    python scripts/run_audit_pipeline.py
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from adora2a_audit.benchmark import (
    analyze_affinity_correlation,
    audit_cognate_redocking_multiseed,
    evaluate_2d_scaffold_vs_random_split,
    plot_affinity_vs_scores,
    plot_metrics_contrast,
    plot_roc_and_early_enrichment,
    prepare_4eiy_hydrated_receptor,
    prepare_4eiy_receptor_and_reference_ligand,
    score_2d_ecfp4_similarity,
    score_3d_ligand_shape,
    score_autodock_vina,
)
from adora2a_audit.curation import audit_reference_ligand_alerts, curate_chembl_activities
from adora2a_audit.decoys import (
    DESCRIPTOR_COLS,
    generate_property_matched_decoys,
    run_1d_bias_audit,
    select_benchmark_actives,
)
from adora2a_audit.metrics import bootstrap_screening_ci, evaluate_screening
from adora2a_audit.structures import (
    build_structure_audit_table,
    inspect_pocket_coordinates,
    inspect_pocket_waters,
    plot_structure_audit,
)

ROOT = Path(__file__).resolve().parents[1]
CACHE_DIR = ROOT / "data" / "cache"
PROC_DIR = ROOT / "data" / "processed"
FIG_DIR = ROOT / "reports" / "figures"
VINA_EXE = CACHE_DIR / "bin" / "vina.exe"


def plot_decoy_bias_control(
    actives_df: pd.DataFrame,
    matched_df: pd.DataFrame,
    unmatched_df: pd.DataFrame,
    bias_res,
    out_path: Path,
) -> None:
    """Visualize 1D physicochemical distributions (Matched vs Unmatched) and go/no-go gate."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(14.0, 4.2), dpi=180)
    cols_to_show = [
        ("mw", "Molecular Weight (Da)"),
        ("clogp", "cLogP (Wildman-Crippen)"),
        ("tpsa", "TPSA (Å²) — Typical GPCR Bias"),
    ]
    for ax, (col, title) in zip(axes[:3], cols_to_show):
        bins = np.linspace(
            min(actives_df[col].min(), unmatched_df[col].min()),
            max(actives_df[col].max(), unmatched_df[col].max()),
            18,
        )
        ax.hist(
            unmatched_df[col],
            bins=bins,
            density=True,
            alpha=0.35,
            color="#ef4444",
            label="Random decoys (unmatched)",
        )
        ax.hist(
            matched_df[col],
            bins=bins,
            density=True,
            alpha=0.45,
            color="#3b82f6",
            label="7D-matched decoys",
        )
        ax.hist(
            actives_df[col],
            bins=bins,
            density=True,
            histtype="step",
            lw=2.2,
            color="#111827",
            label="ADORA2A actives",
        )
        ax.set_title(title, fontsize=9.5, fontweight="bold")
        ax.set_ylabel("Density", fontsize=8.5)
        ax.grid(True, alpha=0.2)
        if col == "mw":
            ax.legend(loc="upper right", fontsize=7.0)

    ax4 = axes[3]
    labels = ["LogReg\n(Random)", "RF\n(Random)", "LogReg\n(7D Match)", "RF\n(7D Match)"]
    aucs = [
        bias_res.unmatched_logistic_auc,
        bias_res.unmatched_rf_auc,
        bias_res.logistic_metrics.roc_auc,
        bias_res.rf_metrics.roc_auc,
    ]
    bar_cols = ["#fca5a5", "#ef4444", "#86efac", "#22c55e"]
    bars = ax4.bar(labels, aucs, color=bar_cols, edgecolor="#334155", width=0.58)
    ax4.axhline(0.60, color="#dc2626", linestyle="--", lw=1.2, label="NO-GO cutoff (0.60)")
    ax4.axhline(0.50, color="#64748b", linestyle=":", lw=1.0, label="Ideal random (0.50)")
    for b, val in zip(bars, aucs, strict=False):
        ax4.text(
            b.get_x() + b.get_width() / 2,
            val + 0.02,
            f"{val:.2f}",
            ha="center",
            fontsize=8.5,
            fontweight="bold",
        )
    ax4.set_ylim(0, 0.98)
    ax4.set_ylabel("ROC-AUC (5-fold CV, 1D props only)", fontsize=8.5)
    ax4.set_title("1D Leakage Go/No-Go Audit", fontsize=9.5, fontweight="bold")
    ax4.tick_params(axis="x", labelsize=7.5)
    ax4.legend(loc="upper right", fontsize=7.2)
    ax4.grid(True, axis="y", alpha=0.2)

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    t_start = time.perf_counter()
    PROC_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("=== PHASE 1: GPCRdb + RCSB PDB Structural & Pocket Water Audit ===")
    struct_df = build_structure_audit_table(CACHE_DIR)
    struct_df.to_csv(PROC_DIR / "adora2a_structures_audit.csv", index=False)
    plot_structure_audit(struct_df, FIG_DIR / "structure_audit_resolution_state.png")

    p_4eiy = inspect_pocket_coordinates(
        CACHE_DIR / "structures" / "pdb_00004eiy.cif", pdb_code="4EIY", ligand_resname="ZMA"
    )
    p_5iu4 = inspect_pocket_coordinates(
        CACHE_DIR / "structures" / "pdb_00005iu4.cif", pdb_code="5IU4", ligand_resname="ZMA"
    )
    p_9ee8 = inspect_pocket_coordinates(
        CACHE_DIR / "structures" / "pdb_00009ee8.cif", pdb_code="9EE8", ligand_resname="ADN"
    )

    pdb_4eiy_path = CACHE_DIR / "structures" / "pdb_00004eiy.pdb"
    pocket_waters = inspect_pocket_waters(pdb_4eiy_path)
    waters_df = pd.DataFrame(
        [
            {
                **asdict(w),
                "xyz": f"({w.xyz[0]:.2f}, {w.xyz[1]:.2f}, {w.xyz[2]:.2f})",
                "receptor_hbonds": "; ".join(w.receptor_hbonds),
                "water_hbonds": "; ".join(str(x) for x in w.water_hbonds),
            }
            for w in pocket_waters
        ]
    )
    waters_df.to_csv(PROC_DIR / "4eiy_pocket_waters_audit.csv", index=False)
    retained_waters = [w for w in pocket_waters if w.retained]

    print(f"  Total ADORA2A structures audited: {len(struct_df)}")
    print(f"  4EIY WT pocket verified: {p_4eiy.is_wt_orthosteric_pocket} | ZMA box center: {p_4eiy.ligand_centroid_xyz}")
    print(
        f"  4EIY crystallographic waters within 5.0 A: {len(pocket_waters)} -> "
        f"Retained ({len(retained_waters)}): {[w.seqid for w in retained_waters]}"
    )

    print("\n=== PHASE 2: ChEMBL251 Bioactivity Curation & Decoy Bias Audit ===")
    raw_prov_df, dedup_ki_df, waterfall = curate_chembl_activities(CACHE_DIR)
    raw_prov_df.to_csv(PROC_DIR / "chembl251_curated_measurements.csv", index=False)
    dedup_ki_df.to_csv(PROC_DIR / "chembl251_deduplicated_ki.csv", index=False)

    ref_alerts_df = audit_reference_ligand_alerts()
    ref_alerts_df.to_csv(PROC_DIR / "reference_ligands_alerts_audit.csv", index=False)

    actives_df = select_benchmark_actives(dedup_ki_df, n_actives=25, min_pki=7.8, seed=42)
    all_adora2a_inchikeys = set(raw_prov_df["inchikey"].unique().tolist())
    matched_decoys_df, unmatched_decoys_df = generate_property_matched_decoys(
        actives_df,
        CACHE_DIR,
        excluded_inchikeys=all_adora2a_inchikeys,
        decoys_per_active=39,
        max_tanimoto=0.35,
        seed=42,
    )

    bias_res, _, oof_1d_rf = run_1d_bias_audit(
        actives_df, matched_decoys_df, unmatched_decoys_df, max_auc_threshold=0.60, seed=42
    )
    print(
        f"  1D Bias Gate: Unmatched AUC (LR={bias_res.unmatched_logistic_auc:.3f}, RF={bias_res.unmatched_rf_auc:.3f}) "
        f"-> Matched AUC (LR={bias_res.logistic_metrics.roc_auc:.3f}, RF={bias_res.rf_metrics.roc_auc:.3f}) "
        f"| PASSED GO/NO-GO: {bias_res.passed_go_nogo}"
    )
    if not bias_res.passed_go_nogo:
        raise RuntimeError("Benchmark failed 1D physicochemical bias gate (ROC-AUC > 0.60)!")

    plot_decoy_bias_control(
        actives_df,
        matched_decoys_df,
        unmatched_decoys_df,
        bias_res,
        FIG_DIR / "decoy_property_matching_audit.png",
    )

    print("\n=== PHASE 3: Clean A/B Benchmark (4EIY Dry vs. 4EIY Hydrated) & 10-Seed Redocking ===")
    vina_prep_dir = CACHE_DIR / "vina_prep"
    rec_dry_pdbqt, zma_3d_h, box_center = prepare_4eiy_receptor_and_reference_ligand(
        pdb_4eiy_path, vina_prep_dir
    )
    rec_hyd_pdbqt = prepare_4eiy_hydrated_receptor(
        rec_dry_pdbqt, retained_waters, vina_prep_dir
    )

    print("  Running 10-seed cognate ZM241385 redocking on 4EIY_dry and 4EIY_hydrated...")
    redock_dry = audit_cognate_redocking_multiseed(
        VINA_EXE, rec_dry_pdbqt, box_center, pdb_4eiy_path, receptor_arm="4EIY_dry"
    )
    redock_hyd = audit_cognate_redocking_multiseed(
        VINA_EXE, rec_hyd_pdbqt, box_center, pdb_4eiy_path, receptor_arm="4EIY_hydrated"
    )
    print(
        f"  [10-Seed Redocking · DRY] Pose1 RMSD={redock_dry.pose1_rmsd_mean:.2f}±{redock_dry.pose1_rmsd_std:.2f} A "
        f"({redock_dry.pose1_lt_2a_count}/10 < 2.0 A) | Best RMSD={redock_dry.best_rmsd_mean:.2f}±{redock_dry.best_rmsd_std:.2f} A "
        f"({redock_dry.any_pose_lt_2a_count}/10 < 2.0 A) | ΔE={redock_dry.delta_e_native_minus_nonnative_mean:+.3f}±{redock_dry.delta_e_native_minus_nonnative_std:.3f} kcal/mol"
    )
    print(
        f"  [10-Seed Redocking · HYDRATED] Pose1 RMSD={redock_hyd.pose1_rmsd_mean:.2f}±{redock_hyd.pose1_rmsd_std:.2f} A "
        f"({redock_hyd.pose1_lt_2a_count}/10 < 2.0 A) | Best RMSD={redock_hyd.best_rmsd_mean:.2f}±{redock_hyd.best_rmsd_std:.2f} A "
        f"({redock_hyd.any_pose_lt_2a_count}/10 < 2.0 A) | ΔE={redock_hyd.delta_e_native_minus_nonnative_mean:+.3f}±{redock_hyd.delta_e_native_minus_nonnative_std:.3f} kcal/mol"
    )

    # Assemble 1,000-compound retrospective set (25 actives + 975 property-matched decoys)
    common_cols = [
        "inchikey",
        "canonical_smiles",
        "murcko_smiles",
        "primary_chembl_id",
        *DESCRIPTOR_COLS,
        "is_active",
    ]
    act_part = actives_df[common_cols + ["pki_median", "reference_name"]].copy()
    dec_part = matched_decoys_df[common_cols].copy()
    dec_part["pki_median"] = np.nan
    dec_part["reference_name"] = ""
    bench_df = pd.concat([act_part, dec_part], ignore_index=True)
    bench_df["score_1d_rf"] = oof_1d_rf
    bench_df["score_2d_zma"] = score_2d_ecfp4_similarity(bench_df["canonical_smiles"].tolist())

    cached_scores_path = PROC_DIR / "benchmark_1000_scores.csv"
    prev_bench = pd.read_csv(cached_scores_path) if cached_scores_path.exists() else None
    same_bench = prev_bench is not None and list(prev_bench["inchikey"]) == list(bench_df["inchikey"])

    if same_bench and {"score_3d_usrcat", "score_3d_shape"}.issubset(prev_bench.columns):
        print("  Reusing cached 3D USRCAT & Shape Tanimoto scores from benchmark_1000_scores.csv")
        bench_df["score_3d_usrcat"] = prev_bench["score_3d_usrcat"].to_numpy(dtype=float)
        bench_df["score_3d_shape"] = prev_bench["score_3d_shape"].to_numpy(dtype=float)
    else:
        print("  Computing 3D conformers + USRCAT & CrippenO3A Shape Tanimoto (1,000 compounds)...")
        usr_scores, shape_combo = score_3d_ligand_shape(bench_df["canonical_smiles"].tolist(), zma_3d_h)
        bench_df["score_3d_usrcat"] = usr_scores
        bench_df["score_3d_shape"] = shape_combo

    if same_bench and "vina_kcal" in prev_bench.columns:
        print("  Reusing cached AutoDock Vina (4EIY dry) scores from benchmark_1000_scores.csv")
        bench_df["vina_kcal"] = prev_bench["vina_kcal"].to_numpy(dtype=float)
    else:
        print("  Running AutoDock Vina v1.2.7 (4EIY dry) across 1,000 benchmark molecules...")
        bench_df["vina_kcal"] = score_autodock_vina(
            bench_df["canonical_smiles"].tolist(), VINA_EXE, rec_dry_pdbqt, box_center, tag="4EIY dry"
        )
    bench_df["vina_neg_kcal"] = -bench_df["vina_kcal"]

    if same_bench and "vina_hyd_kcal" in prev_bench.columns:
        print("  Reusing cached AutoDock Vina (4EIY hydrated) scores from benchmark_1000_scores.csv")
        bench_df["vina_hyd_kcal"] = prev_bench["vina_hyd_kcal"].to_numpy(dtype=float)
    else:
        print("  Running AutoDock Vina v1.2.7 (4EIY hydrated) across 1,000 benchmark molecules...")
        bench_df["vina_hyd_kcal"] = score_autodock_vina(
            bench_df["canonical_smiles"].tolist(), VINA_EXE, rec_hyd_pdbqt, box_center, tag="4EIY hydrated"
        )
    bench_df["vina_hyd_neg_kcal"] = -bench_df["vina_hyd_kcal"]
    bench_df.to_csv(cached_scores_path, index=False)

    # Wide-dynamic-range experimental pKi panel (n = 60 measured binders)
    print("  Evaluating wide-range experimental pKi panel (n = 60 measured Ki binders)...")
    non_ago_ki = dedup_ki_df[
        (dedup_ki_df["modality"] != "Agonist")
        & (~dedup_ki_df["inter_assay_discrepancy"])
        & (dedup_ki_df["mw"].between(180.0, 540.0))
    ].copy()
    already_ik = set(actives_df["inchikey"])
    extra_refs = non_ago_ki[
        (non_ago_ki["reference_name"] != "") & (~non_ago_ki["inchikey"].isin(already_ik))
    ]
    low_mid_pool = non_ago_ki[
        (~non_ago_ki["inchikey"].isin(already_ik))
        & (~non_ago_ki["inchikey"].isin(set(extra_refs["inchikey"])))
        & (non_ago_ki["pki_median"] < 7.8)
    ].sort_values("pki_median")
    n_extra_target = max(0, 35 - len(extra_refs))
    lin_idx = np.linspace(0, len(low_mid_pool) - 1, n_extra_target, dtype=int)
    extra_panel = pd.concat([extra_refs, low_mid_pool.iloc[lin_idx]], ignore_index=True)

    extra_cache = PROC_DIR / "affinity_panel_60_scores.csv"
    prev_ep = pd.read_csv(extra_cache) if extra_cache.exists() else None
    panel_cols = [
        "inchikey",
        "canonical_smiles",
        "primary_chembl_id",
        "reference_name",
        "pki_median",
        "vina_kcal",
        "vina_neg_kcal",
        "vina_hyd_kcal",
        "vina_hyd_neg_kcal",
        "score_2d_zma",
        "score_3d_shape",
    ]

    if prev_ep is not None and len(prev_ep) == len(actives_df) + len(extra_panel):
        if "vina_hyd_kcal" in prev_ep.columns:
            affinity_panel = prev_ep
        else:
            # Reuse dry Vina + 2D + 3D for the 35 extra compounds, dock only hydrated Vina on the 35 extra!
            extra_prev = prev_ep.iloc[len(actives_df) :].copy().reset_index(drop=True)
            extra_prev["vina_hyd_kcal"] = score_autodock_vina(
                extra_prev["canonical_smiles"].tolist(),
                VINA_EXE,
                rec_hyd_pdbqt,
                box_center,
                tag="4EIY hydrated (panel-35)",
            )
            extra_prev["vina_hyd_neg_kcal"] = -extra_prev["vina_hyd_kcal"]
            affinity_panel = pd.concat(
                [bench_df[bench_df["is_active"] == 1][panel_cols], extra_prev[panel_cols]],
                ignore_index=True,
            )
            affinity_panel.to_csv(extra_cache, index=False)
    else:
        extra_panel["score_2d_zma"] = score_2d_ecfp4_similarity(extra_panel["canonical_smiles"].tolist())
        _, extra_shape = score_3d_ligand_shape(extra_panel["canonical_smiles"].tolist(), zma_3d_h)
        extra_panel["score_3d_shape"] = extra_shape
        extra_panel["vina_kcal"] = score_autodock_vina(
            extra_panel["canonical_smiles"].tolist(), VINA_EXE, rec_dry_pdbqt, box_center, tag="4EIY dry (panel-35)"
        )
        extra_panel["vina_neg_kcal"] = -extra_panel["vina_kcal"]
        extra_panel["vina_hyd_kcal"] = score_autodock_vina(
            extra_panel["canonical_smiles"].tolist(), VINA_EXE, rec_hyd_pdbqt, box_center, tag="4EIY hydrated (panel-35)"
        )
        extra_panel["vina_hyd_neg_kcal"] = -extra_panel["vina_hyd_kcal"]
        affinity_panel = pd.concat(
            [bench_df[bench_df["is_active"] == 1][panel_cols], extra_panel[panel_cols]],
            ignore_index=True,
        )
        affinity_panel.to_csv(extra_cache, index=False)

    # Murcko Scaffold Split vs Stratified Random Split on experimental CHEMBL251 Ki data
    exp_split_df = non_ago_ki[
        (non_ago_ki["pki_median"] >= 9.0) | (non_ago_ki["pki_median"] <= 6.5)
    ].copy()
    exp_y = (exp_split_df["pki_median"] >= 9.0).astype(int).to_numpy()
    ml_rand_metrics, ml_scaf_metrics = evaluate_2d_scaffold_vs_random_split(
        exp_split_df["canonical_smiles"].tolist(),
        exp_split_df["murcko_smiles"].tolist(),
        exp_y,
        seed=42,
        label_prefix="2D ECFP4 ML experimental",
    )

    # Evaluate all 5 benchmark series + 95% bootstrap CIs
    y_bench = bench_df["is_active"].to_numpy(dtype=int)
    method_scores = {
        "1D Physicochemical (Control)": bench_df["score_1d_rf"].to_numpy(dtype=float),
        "2D ECFP4 (vs ZM241385)": bench_df["score_2d_zma"].to_numpy(dtype=float),
        "3D Shape + USRCAT (vs ZMA)": bench_df["score_3d_shape"].to_numpy(dtype=float),
        "3D Docking (Vina · 4EIY dry)": bench_df["vina_neg_kcal"].to_numpy(dtype=float),
        "3D Docking (Vina · 4EIY hydrated)": bench_df["vina_hyd_neg_kcal"].to_numpy(dtype=float),
    }

    non_analogue_mask = (bench_df["is_active"] == 0) | (bench_df["score_2d_zma"] < 0.35)
    y_hop = bench_df.loc[non_analogue_mask, "is_active"].to_numpy(dtype=int)
    hop_2d = evaluate_screening(
        y_hop,
        bench_df.loc[non_analogue_mask, "score_2d_zma"].to_numpy(dtype=float),
        method="2D ECFP4 (no ZMA analogues, ECFP4<0.35)",
    )
    hop_vina_dry = evaluate_screening(
        y_hop,
        bench_df.loc[non_analogue_mask, "vina_neg_kcal"].to_numpy(dtype=float),
        method="AutoDock Vina · 4EIY dry (no ZMA analogues)",
    )
    hop_vina_hyd = evaluate_screening(
        y_hop,
        bench_df.loc[non_analogue_mask, "vina_hyd_neg_kcal"].to_numpy(dtype=float),
        method="AutoDock Vina · 4EIY hydrated (no ZMA analogues)",
    )

    metrics_list = []
    ci_map = {}
    for name, sc in method_scores.items():
        m = evaluate_screening(y_bench, sc, method=name)
        ci = bootstrap_screening_ci(y_bench, sc, n_bootstrap=500, seed=42)
        metrics_list.append(m)
        ci_map[name] = ci
        print(
            f"  [{name}] EF1%={m.ef1:.2f}x (CI {ci['ef1'][0]:.1f}-{ci['ef1'][1]:.1f}), "
            f"EF5%={m.ef5:.2f}x (CI {ci['ef5'][0]:.1f}-{ci['ef5'][1]:.1f}), "
            f"BEDROC20={m.bedroc_20:.3f}, ROC-AUC={m.roc_auc:.3f}, PR-AUC={m.pr_auc:.3f}, Acc@1%={m.accuracy_top1pct:.3f}"
        )

    print(
        f"  [Scaffold-hop subset, n_actives={int(y_hop.sum())}] "
        f"2D ECFP4: EF1%={hop_2d.ef1:.2f}x, AUC={hop_2d.roc_auc:.3f} | "
        f"Vina dry: EF1%={hop_vina_dry.ef1:.2f}x, AUC={hop_vina_dry.roc_auc:.3f} | "
        f"Vina hydrated: EF1%={hop_vina_hyd.ef1:.2f}x, AUC={hop_vina_hyd.roc_auc:.3f}"
    )

    corr_potent = analyze_affinity_correlation(bench_df[bench_df["is_active"] == 1])
    corr_wide = analyze_affinity_correlation(affinity_panel)
    print("  Affinity correlation (wide range n=60):", json.dumps(corr_wide, indent=2))

    plot_roc_and_early_enrichment(
        y_bench, method_scores, FIG_DIR / "benchmark_roc_and_enrichment.png"
    )
    plot_metrics_contrast(
        metrics_list, ci_map, FIG_DIR / "metrics_accuracy_vs_enrichment.png"
    )
    scatter_df = pd.concat(
        [
            affinity_panel.assign(is_active=1),
            bench_df[bench_df["is_active"] == 0][panel_cols].assign(is_active=0),
        ],
        ignore_index=True,
    )
    plot_affinity_vs_scores(
        scatter_df, corr_wide, FIG_DIR / "affinity_vs_score_scatter.png"
    )

    summary_payload = {
        "n_structures_audited": len(struct_df),
        "structures_by_state": struct_df["state"].value_counts().to_dict(),
        "structures_na_by_state": struct_df.groupby("state")["has_allosteric_na"].sum().astype(int).to_dict(),
        "pocket_checks": {
            "4EIY_wt": p_4eiy.is_wt_orthosteric_pocket,
            "4EIY_center": p_4eiy.ligand_centroid_xyz,
            "5IU4_mutations": {str(k): list(v) for k, v in p_5iu4.mutated_pocket_residues.items()},
            "9EE8_wt": p_9ee8.is_wt_orthosteric_pocket,
        },
        "pocket_waters_audit": {
            "n_within_5a": len(pocket_waters),
            "n_retained": len(retained_waters),
            "retained_seqids": [w.seqid for w in retained_waters],
            "excluded_seqids": [w.seqid for w in pocket_waters if not w.retained],
            "criterion": "min_dist_to_ZMA <= 4.5 A, B-factor < 30.0 A^2, occupancy == 1.0, >=1 H-bond <= 3.5 A",
        },
        "curation_waterfall": waterfall,
        "bias_audit": {
            "unmatched_logistic_auc": round(bias_res.unmatched_logistic_auc, 4),
            "unmatched_rf_auc": round(bias_res.unmatched_rf_auc, 4),
            "matched_logistic_auc": round(bias_res.logistic_metrics.roc_auc, 4),
            "matched_rf_auc": round(bias_res.rf_metrics.roc_auc, 4),
            "passed_go_nogo": bias_res.passed_go_nogo,
        },
        "cognate_redocking_10seeds": {
            "dry": redock_dry.to_dict(),
            "hydrated": redock_hyd.to_dict(),
        },
        "benchmark_metrics": [m.to_dict() for m in metrics_list],
        "benchmark_bootstrap_ci": ci_map,
        "scaffold_hop_no_zma_analogues": {
            "n_actives_remaining": int(y_hop.sum()),
            "2d_ecfp4": hop_2d.to_dict(),
            "autodock_vina_dry": hop_vina_dry.to_dict(),
            "autodock_vina_hydrated": hop_vina_hyd.to_dict(),
        },
        "murcko_split_experimental_chembl": {
            "n_compounds": len(exp_split_df),
            "random_split": ml_rand_metrics.to_dict(),
            "murcko_scaffold_split": ml_scaf_metrics.to_dict(),
        },
        "affinity_correlation_potent_n25": corr_potent,
        "affinity_correlation_wide_n60": corr_wide,
        "runtime_seconds": round(time.perf_counter() - t_start, 1),
    }
    (PROC_DIR / "audit_summary.json").write_text(
        json.dumps(summary_payload, indent=2), encoding="utf-8"
    )
    print(f"\nPipeline v2 completed cleanly in {summary_payload['runtime_seconds']}s!")


if __name__ == "__main__":
    main()
