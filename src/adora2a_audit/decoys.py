"""Property-matched decoy selection and 1D physicochemical bias audit (go/no-go gate).

Selects real drug-like ChEMBL molecules matched jointly in standardized 7D
physicochemical space (MW, cLogP, TPSA, HBD, HBA, RotatableBonds, FormalCharge)
while enforcing topological dissimilarity (Morgan/ECFP4 Tanimoto < 0.35 to ALL actives)
and excluding any molecule ever tested on ADORA2A (`CHEMBL251`).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from adora2a_audit.curation import standardize_smiles
from adora2a_audit.metrics import ScreeningMetrics, evaluate_screening

DESCRIPTOR_COLS = [
    "mw",
    "clogp",
    "tpsa",
    "hbd",
    "hba",
    "rot_bonds",
    "formal_charge",
]

_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


@dataclass(frozen=True)
class BiasAuditResult:
    """Go/no-go audit result for 1D physicochemical benchmark leakage."""

    logistic_metrics: ScreeningMetrics
    rf_metrics: ScreeningMetrics
    unmatched_logistic_auc: float
    unmatched_rf_auc: float
    passed_go_nogo: bool
    max_auc_threshold: float = 0.60


def compute_1d_descriptors_and_fp(smiles: str):
    """Compute the 7 physicochemical descriptors and Morgan radius-2 (ECFP4) fingerprint."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    props = {
        "mw": float(Descriptors.MolWt(mol)),
        "clogp": float(Crippen.MolLogP(mol)),
        "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "rot_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "formal_charge": int(Chem.GetFormalCharge(mol)),
    }
    fp = _MORGAN_GEN.GetFingerprint(mol)
    return props, fp


def select_benchmark_actives(
    dedup_ki_df: pd.DataFrame,
    n_actives: int = 50,
    min_pki: float = 7.5,
    mw_range: tuple[float, float] = (220.0, 520.0),
    seed: int = 42,
) -> pd.DataFrame:
    """Select a chemically diverse, high-confidence antagonist/non-agonist active set.

    Prioritizes:
    1. Canonical co-crystallized / clinical ADORA2A antagonists present in curated Ki_B
       (e.g. ZM241385, Istradefylline, Preladenant, Vipadenant, Tozadenant, Imaradenant, Theophylline/Caffeine if potent, etc.)
    2. Murcko-scaffold-diverse potent actives (`pKi >= min_pki`, no >1.0 log assay discrepancy,
       non-agonist modality, drug-like MW range) so the benchmark tests genuine scaffold diversity
       rather than 50 methyl/chloro analogues of a single congeneric series.
    """
    elig = dedup_ki_df[
        (dedup_ki_df["modality"] != "Agonist")
        & (~dedup_ki_df["inter_assay_discrepancy"])
        & (dedup_ki_df["pki_median"] >= min_pki)
        & (dedup_ki_df["mw"].between(mw_range[0], mw_range[1]))
    ].copy()

    # Always include reference antagonists that meet eligibility
    ref_Part = elig[elig["reference_name"] != ""].copy()
    rest = elig[elig["reference_name"] == ""].sort_values("pki_median", ascending=False).copy()

    selected_indices: list[int] = list(ref_Part.index)
    seen_scaffolds: set[str] = set(ref_Part["murcko_smiles"].tolist())

    # Pass 1: one representative (most potent) per distinct Murcko scaffold
    for idx, row in rest.iterrows():
        if len(selected_indices) >= n_actives:
            break
        scaf = row["murcko_smiles"]
        if scaf not in seen_scaffolds:
            selected_indices.append(idx)
            seen_scaffolds.add(scaf)

    # Pass 2: if still under n_actives, fill deterministically from remaining potent actives
    if len(selected_indices) < n_actives:
        rng = np.random.default_rng(seed)
        remaining = [idx for idx in rest.index if idx not in set(selected_indices)]
        extra = rng.choice(remaining, size=min(n_actives - len(selected_indices), len(remaining)), replace=False)
        selected_indices.extend(extra.tolist())

    out = elig.loc[selected_indices].copy().reset_index(drop=True)
    desc_rows = []
    for smi in out["canonical_smiles"]:
        props, _ = compute_1d_descriptors_and_fp(smi)
        desc_rows.append(props)
    desc_df = pd.DataFrame(desc_rows)
    for col in DESCRIPTOR_COLS:
        out[col] = desc_df[col]
    out["is_active"] = 1
    return out


def generate_property_matched_decoys(
    actives_df: pd.DataFrame,
    cache_dir: str | Path,
    excluded_inchikeys: set[str],
    decoys_per_active: int = 19,
    max_tanimoto: float = 0.35,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select `decoys_per_active` unique property-matched decoys per active from ChEMBL pool.

    Also returns an equal-sized unmatched random decoy sample from the same pool to
    demonstrate how severe 1D physicochemical bias is when property matching is omitted.
    """
    cache_path = Path(cache_dir)
    pool_raw = json.loads((cache_path / "chembl_decoy_pool.json").read_text(encoding="utf-8"))
    excluded_conn = {ik[:14] for ik in excluded_inchikeys}

    active_fps = []
    for smi in actives_df["canonical_smiles"]:
        _, fp = compute_1d_descriptors_and_fp(smi)
        active_fps.append(fp)

    # Standardize and filter candidate pool
    cand_records = []
    seen_conn = set(excluded_conn)
    for item in pool_raw:
        std = standardize_smiles(item.get("canonical_smiles", ""))
        if std is None:
            continue
        conn = std.inchikey[:14]
        if conn in seen_conn:
            continue
        props, fp = compute_1d_descriptors_and_fp(std.canonical_smiles)
        if props is None or fp is None:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp, active_fps)
        max_sim = float(max(sims)) if sims else 0.0
        if max_sim >= max_tanimoto:
            continue
        seen_conn.add(conn)
        cand_records.append(
            {
                "inchikey": std.inchikey,
                "canonical_smiles": std.canonical_smiles,
                "murcko_smiles": std.murcko_smiles,
                "primary_chembl_id": item.get("molecule_chembl_id", ""),
                "max_tanimoto_to_actives": round(max_sim, 4),
                **props,
            }
        )

    cand_df = pd.DataFrame(cand_records)
    n_needed = len(actives_df) * decoys_per_active
    if len(cand_df) < n_needed:
        raise RuntimeError(
            f"Candidate pool ({len(cand_df)}) smaller than requested decoys ({n_needed})."
        )

    # Standardize descriptor space using active-set mean and std (with floor to prevent zero-div)
    act_mat = actives_df[DESCRIPTOR_COLS].to_numpy(dtype=float)
    cand_mat = cand_df[DESCRIPTOR_COLS].to_numpy(dtype=float)
    mu = act_mat.mean(axis=0)
    sigma = np.maximum(act_mat.std(axis=0), np.array([15.0, 0.4, 10.0, 0.5, 0.8, 0.8, 0.25]))

    # Emphasize TPSA, MW, cLogP, HBD, HBA equally so nonlinear 1D trees cannot exploit TPSA/HBA gaps
    weights = np.array([1.2, 1.2, 1.5, 1.2, 1.3, 1.0, 2.0])
    act_z = ((act_mat - mu) / sigma) * weights
    cand_z = ((cand_mat - mu) / sigma) * weights

    # Greedy round-robin assignment (1 decoy per active per round across `decoys_per_active` rounds)
    # ensures every active gets fair first-choice picks without order starvation.
    used = np.zeros(len(cand_df), dtype=bool)
    matched_rows = []
    rng = np.random.default_rng(seed)
    active_order = np.arange(len(actives_df))

    for r in range(decoys_per_active):
        shuffled_actives = rng.permutation(active_order)
        for a_idx in shuffled_actives:
            diffs = cand_z - act_z[a_idx]
            dists = np.einsum("ij,ij->i", diffs, diffs)
            dists[used] = np.inf
            best_c = int(np.argmin(dists))
            used[best_c] = True
            row_dict = cand_df.iloc[best_c].to_dict()
            row_dict["matched_active_inchikey"] = actives_df.iloc[a_idx]["inchikey"]
            row_dict["prop_z_distance"] = round(float(np.sqrt(dists[best_c])), 4)
            row_dict["is_active"] = 0
            matched_rows.append(row_dict)

    matched_df = pd.DataFrame(matched_rows).reset_index(drop=True)

    # Also draw an unmatched random decoy sample of the same size from cand_df
    rand_idx = rng.choice(len(cand_df), size=n_needed, replace=False)
    unmatched_df = cand_df.iloc[rand_idx].copy().reset_index(drop=True)
    unmatched_df["is_active"] = 0

    return matched_df, unmatched_df


def run_1d_bias_audit(
    actives_df: pd.DataFrame,
    matched_decoys_df: pd.DataFrame,
    unmatched_decoys_df: pd.DataFrame,
    max_auc_threshold: float = 0.60,
    seed: int = 42,
) -> tuple[BiasAuditResult, np.ndarray, np.ndarray]:
    """Evaluate 5-fold cross-validated 1D physicochemical classifiers on matched vs unmatched sets.

    Returns:
        audit_result: Summary with go/no-go flag (`passed_go_nogo` requires BOTH LR and RF < threshold)
        y_true_matched: Combined labels `[1..1, 0..0]`
        oof_probs_rf: Out-of-fold 1D RandomForest probabilities for downstream benchmark comparison
    """
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    lr = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000, random_state=seed))
    rf = RandomForestClassifier(
        n_estimators=200,
        max_depth=4,
        min_samples_leaf=8,
        random_state=seed,
    )

    def _eval_pair(dec_df: pd.DataFrame, tag: str):
        comb = pd.concat(
            [actives_df[DESCRIPTOR_COLS + ["is_active"]], dec_df[DESCRIPTOR_COLS + ["is_active"]]],
            ignore_index=True,
        )
        x = comb[DESCRIPTOR_COLS].to_numpy(dtype=float)
        y = comb["is_active"].to_numpy(dtype=int)
        p_lr = cross_val_predict(lr, x, y, cv=cv, method="predict_proba")[:, 1]
        p_rf = cross_val_predict(rf, x, y, cv=cv, method="predict_proba")[:, 1]
        m_lr = evaluate_screening(y, p_lr, method=f"1D Logistic ({tag})")
        m_rf = evaluate_screening(y, p_rf, method=f"1D RandomForest ({tag})")
        return y, p_lr, p_rf, m_lr, m_rf

    y_m, _, p_rf_m, m_lr_m, m_rf_m = _eval_pair(matched_decoys_df, "matched")
    _, _, _, m_lr_u, m_rf_u = _eval_pair(unmatched_decoys_df, "unmatched")

    passed = (m_lr_m.roc_auc <= max_auc_threshold) and (m_rf_m.roc_auc <= max_auc_threshold)
    result = BiasAuditResult(
        logistic_metrics=m_lr_m,
        rf_metrics=m_rf_m,
        unmatched_logistic_auc=m_lr_u.roc_auc,
        unmatched_rf_auc=m_rf_u.roc_auc,
        passed_go_nogo=passed,
        max_auc_threshold=max_auc_threshold,
    )
    return result, y_m, p_rf_m
