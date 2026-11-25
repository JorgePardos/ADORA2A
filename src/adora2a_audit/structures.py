"""GPCRdb + RCSB PDB structural audit for human Adenosine A2A receptor (`ADORA2A`).

Classifies all experimentally determined structures by conformational state
(`Inactive`, `Intermediate`, `Active`), experimental method, resolution,
co-crystallized ligand modality, intracellular fusion partners (`BRIL`, `T4L`,
`Rubredoxin`), signaling complexes (`Gs` / `mini-Gs` / `Nb35` / `Fab`),
engineered receptor mutations (`StaR2` vs. wild-type 7TM pocket), and presence
of the conserved Class A allosteric Na+ ion at Asp52^{2.50}.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import gemmi
import httpx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GPCRDB_URL = "https://gpcrdb.org/services/structure/protein/aa2ar_human/"

# Key human ADORA2A residues (UniProt P29274 numbering) that line the orthosteric
# and Na+ allosteric pockets and must remain wild-type for unbiased docking.
EXPECTED_WT_POCKET: dict[int, str] = {
    52: "ASP",   # 2.50 conserved Na+ allosteric pocket anchor
    88: "THR",   # 3.36 deep ribose / agonist H-bond partner (mutated T88A in StaR2!)
    168: "PHE",  # ECL2 aromatic pi-stacking lid on purine/triazolotriazine core
    169: "GLU",  # ECL2 salt-bridge / exocyclic amine H-bond partner
    246: "TRP",  # 6.48 toggle switch
    250: "HIS",  # 6.52 polar pocket wall
    253: "ASN",  # 6.55 bidentate H-bond anchor conserved across agonists & antagonists
    277: "SER",  # 7.42 ribose H-bond partner (mutated S277A in StaR2!)
    278: "HIS",  # 7.43 polar contact in deep pocket
}


@dataclass(frozen=True)
class PocketInspection:
    """Atomic-level coordinate verification for a candidate reference structure."""

    pdb_code: str
    chain_id: str
    resolved_polymer_residues: int
    observed_pocket_residues: dict[int, str]
    mutated_pocket_residues: dict[int, tuple[str, str]]
    is_wt_orthosteric_pocket: bool
    orthosteric_ligand_id: str
    ligand_centroid_xyz: tuple[float, float, float] | None
    ligand_box_size_xyz: tuple[float, float, float] | None
    has_resolved_na_ion: bool


def fetch_or_load_gpcrdb(cache_dir: str | Path) -> list[dict]:
    """Load cached GPCRdb ADORA2A structures, querying GPCRdb REST API if absent."""
    cache_path = Path(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)
    gpcr_file = cache_path / "gpcrdb_aa2ar_human.json"
    if gpcr_file.exists():
        return json.loads(gpcr_file.read_text(encoding="utf-8"))

    with httpx.Client(timeout=30.0, follow_redirects=True) as client:
        resp = client.get(GPCRDB_URL)
        resp.raise_for_status()
        data = resp.json()
    gpcr_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return data


def build_structure_audit_table(cache_dir: str | Path) -> pd.DataFrame:
    """Merge GPCRdb state/ligand annotations with RCSB PDB construct & mutation metadata."""
    cache_path = Path(cache_dir)
    gpcr_entries = fetch_or_load_gpcrdb(cache_path)
    pdb_meta_path = cache_path / "pdb_metadata_aa2ar.json"
    pdb_meta_map: dict[str, dict] = {}
    if pdb_meta_path.exists():
        raw_meta = json.loads(pdb_meta_path.read_text(encoding="utf-8"))
        pdb_meta_map = {e["rcsb_id"].upper(): e for e in raw_meta.get("entries", [])}

    rows = []
    for g in gpcr_entries:
        pid = g["pdb_code"].upper()
        pm = pdb_meta_map.get(pid, {})
        all_ligs = g.get("ligands") or []
        small_ligs = [
            lig for lig in all_ligs if lig.get("type") == "Small molecule" and lig.get("PDB") != "NA"
        ]
        prot_ligs = [lig for lig in all_ligs if lig.get("type") == "Protein"]

        nonpoly = (pm.get("rcsb_entry_info") or {}).get("nonpolymer_bound_components") or []
        has_na = any(lig.get("PDB") == "NA" for lig in all_ligs) or ("NA" in nonpoly)

        mut_count = 0
        fusions: list[str] = []
        signaling: list[str] = []
        for pe in pm.get("polymer_entities") or []:
            desc = ((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").lower()
            orgs = [
                (o.get("ncbi_scientific_name") or "").lower()
                for o in (pe.get("rcsb_entity_source_organism") or [])
            ]
            if "adenosine" in desc or "a2a" in desc:
                mut_count = max(
                    mut_count, int((pe.get("entity_poly") or {}).get("rcsb_mutation_count") or 0)
                )
                if "b562" in desc or "bril" in desc or "cytochrome" in desc or any("coli" in o for o in orgs):
                    fusions.append("BRIL")
                if "lysozyme" in desc or "endolysin" in desc or any("phage" in o for o in orgs):
                    fusions.append("T4L")
                if "rubredoxin" in desc:
                    fusions.append("Rubredoxin")
            elif any(k in desc for k in ("fab", "antibody", "immunoglobulin", "heavy chain", "light chain")):
                signaling.append("Fab")
            elif any(k in desc for k in ("guanine nucleotide", "mini-gs", "g-alpha", "g protein")):
                signaling.append("Gs/mini-Gs")
            elif "nanobody" in desc or "nb35" in desc:
                signaling.append("Nb35")

        if g.get("signalling_protein") and "Gs/mini-Gs" not in signaling:
            signaling.append("Gs/mini-Gs")

        # Determine primary ligand modality
        funcs = [lig.get("function") or "" for lig in small_ligs]
        if not small_ligs:
            modality = "Antibody/Protein" if prot_ligs else "Apo"
        elif any("Agonist" in f and "Inverse" not in f for f in funcs):
            modality = "Agonist"
        elif any("Inverse" in f for f in funcs):
            modality = "Inverse agonist"
        elif any("Antagonist" in f for f in funcs):
            modality = "Antagonist"
        else:
            modality = "Unknown/Fragment"

        # Classify engineering construct & suitability
        fusion_str = "+".join(sorted(set(fusions))) or "None"
        sig_str = "+".join(sorted(set(signaling))) or "None"
        if mut_count >= 7:
            construct_class = "Thermostabilized StaR2 (pocket T88A/S277A risk)"
        elif mut_count in {4, 5, 6} and fusion_str == "None":
            construct_class = "Thermostabilized GL31/StaR"
        elif "BRIL" in fusions and mut_count <= 3:
            construct_class = "WT-7TM + ICL3-BRIL"
        elif "T4L" in fusions and mut_count <= 3:
            construct_class = "WT-7TM + ICL3-T4L"
        elif "Gs/mini-Gs" in signaling:
            construct_class = "Active Gs/mini-Gs signaling complex"
        else:
            construct_class = "Engineered / Other"

        res = float(g["resolution"]) if g.get("resolution") is not None else np.nan
        pub_date = g.get("publication_date") or ""
        year = int(pub_date[:4]) if len(pub_date) >= 4 and pub_date[:4].isdigit() else np.nan

        rows.append(
            {
                "pdb_code": pid,
                "state": g.get("state") or "Unknown",
                "method": g.get("type") or "Unknown",
                "resolution": res,
                "publication_date": pub_date,
                "year": year,
                "ligand_names": "; ".join(lig.get("name", "") for lig in small_ligs) or "Apo",
                "ligand_pdb_ids": "; ".join(lig.get("PDB", "") for lig in small_ligs if lig.get("PDB")),
                "ligand_modality": modality,
                "ligand_smiles": small_ligs[0].get("SMILES", "") if small_ligs else "",
                "has_allosteric_na": bool(has_na),
                "mutation_count": mut_count,
                "fusion_partner": fusion_str,
                "signaling_partner": sig_str,
                "construct_class": construct_class,
                "unmodeled_residues": (pm.get("rcsb_entry_info") or {}).get(
                    "deposited_unmodeled_polymer_monomer_count"
                ),
                "doi": g.get("publication") or "",
            }
        )

    df = pd.DataFrame(rows).sort_values(["state", "resolution", "pdb_code"]).reset_index(drop=True)
    return df


def inspect_pocket_coordinates(
    structure_path: str | Path,
    pdb_code: str,
    chain_id: str = "A",
    ligand_resname: str | None = None,
    box_padding: float = 8.0,
) -> PocketInspection:
    """Verify wild-type pocket residue identities and extract orthosteric docking box."""
    st = gemmi.read_structure(str(structure_path))
    model = st[0]
    chain = model[chain_id]

    observed: dict[int, str] = {}
    poly_count = 0
    has_na = False
    lig_coords: list[tuple[float, float, float]] = []
    chosen_lig = ligand_resname or ""

    for res in chain:
        seqnum = res.seqid.num
        if res.entity_type == gemmi.EntityType.Polymer:
            poly_count += 1
            if seqnum in EXPECTED_WT_POCKET:
                observed[seqnum] = res.name
        elif res.entity_type == gemmi.EntityType.NonPolymer:
            if res.name == "NA":
                has_na = True
            elif (
                res.name not in {"HOH", "CLR", "OLA", "OLB", "OLC", "PEG", "SO4", "SOG", "GDP"}
                and (ligand_resname is None or res.name == ligand_resname)
            ):
                chosen_lig = res.name
                for atom in res:
                    if atom.element.name != "H":
                        lig_coords.append((atom.pos.x, atom.pos.y, atom.pos.z))

    mutated: dict[int, tuple[str, str]] = {}
    for pos, wt_name in EXPECTED_WT_POCKET.items():
        obs_name = observed.get(pos, "MISSING")
        if obs_name != wt_name:
            mutated[pos] = (wt_name, obs_name)

    centroid = None
    box_size = None
    if lig_coords:
        arr = np.asarray(lig_coords, dtype=float)
        c = arr.mean(axis=0)
        span = arr.max(axis=0) - arr.min(axis=0) + 2.0 * box_padding
        # Enforce sensible Vina box bounds (at least 18 A, at most 24 A per axis)
        span = np.clip(span, 18.0, 24.0)
        centroid = (round(float(c[0]), 3), round(float(c[1]), 3), round(float(c[2]), 3))
        box_size = (round(float(span[0]), 2), round(float(span[1]), 2), round(float(span[2]), 2))

    return PocketInspection(
        pdb_code=pdb_code.upper(),
        chain_id=chain_id,
        resolved_polymer_residues=poly_count,
        observed_pocket_residues=observed,
        mutated_pocket_residues=mutated,
        is_wt_orthosteric_pocket=(len(mutated) == 0),
        orthosteric_ligand_id=chosen_lig,
        ligand_centroid_xyz=centroid,
        ligand_box_size_xyz=box_size,
        has_resolved_na_ion=has_na,
    )


def plot_structure_audit(df: pd.DataFrame, out_path: str | Path) -> None:
    """Create a 2-panel summary plot of all ADORA2A PDB structures (Resolution vs State & Engineering)."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.4), dpi=180)
    palette = {"Inactive": "#1f77b4", "Active": "#d62728", "Intermediate": "#ff7f0e"}
    markers = {True: "o", False: "^"}

    # Panel A: Resolution vs Publication Year by Conformational State & Allosteric Na+
    for state in ["Inactive", "Active", "Intermediate"]:
        sub = df[df["state"] == state]
        for na_flag, label_suffix in [(True, "with allosteric Na+"), (False, "no Na+")]:
            grp = sub[sub["has_allosteric_na"] == na_flag]
            if grp.empty:
                continue
            ax1.scatter(
                grp["year"],
                grp["resolution"],
                c=palette.get(state, "#7f7f7f"),
                marker=markers[na_flag],
                s=52,
                alpha=0.72,
                edgecolors="white",
                linewidths=0.5,
                label=f"{state} ({label_suffix}, n={len(grp)})",
            )

    highlights = {
        "4EIY": ("4EIY (1.80 Å, WT-pocket\nSelected Inactive)", (-52, -32), "#0b4f6c"),
        "5NM4": ("5NM4 (1.70 Å, StaR2: T88A/S277A)", (10, -4), "#555555"),
        "5IU4": ("5IU4 (1.72 Å, StaR2: T88A/S277A)", (10, -24), "#555555"),
        "9EE8": ("9EE8 (2.63 Å, Gs+ADN\nSelected Active)", (-108, -24), "#9b1d20"),
        "4UHR": ("4UHR (2.60 Å, GL31\nNo G protein)", (-58, 18), "#8c564b"),
        "5G53": ("5G53 (3.40 Å, mini-Gs)", (-50, 16), "#9b1d20"),
    }
    for pid, (txt, offset, col) in highlights.items():
        row = df[df["pdb_code"] == pid]
        if not row.empty:
            r = row.iloc[0]
            ax1.scatter([r["year"]], [r["resolution"]], s=110, facecolors="none", edgecolors=col, linewidths=2.0, zorder=6)
            ax1.annotate(
                txt,
                xy=(r["year"], r["resolution"]),
                xytext=offset,
                textcoords="offset points",
                fontsize=7.5,
                fontweight="bold",
                color=col,
                bbox={"boxstyle": "round,pad=0.22", "fc": "white", "ec": col, "alpha": 0.92},
                arrowprops={"arrowstyle": "->", "color": col, "lw": 0.9},
                zorder=7,
            )

    ax1.axhline(2.5, color="#888888", linestyle="--", linewidth=0.9, alpha=0.7)
    ax1.text(2008.3, 2.44, "High-resolution cutoff (2.5 Å)", fontsize=7.5, color="#555555")
    ax1.invert_yaxis()
    ax1.set_xlabel("PDB Release Year", fontsize=10)
    ax1.set_ylabel("Experimental Resolution (Å, higher = better)", fontsize=10)
    ax1.set_title("A. ADORA2A Structural Landscape (n = 188 PDB entries)", fontsize=10.5, fontweight="bold")
    ax1.legend(loc="lower left", fontsize=7.5, framealpha=0.92)
    ax1.grid(True, alpha=0.22)

    # Panel B: Engineering construct classes across conformational states
    summary = (
        df.groupby(["construct_class", "state"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=["Inactive", "Active", "Intermediate"], fill_value=0)
    )
    summary["total"] = summary.sum(axis=1)
    summary = summary.sort_values("total", ascending=True).drop(columns=["total"])

    y_pos = np.arange(len(summary))
    left = np.zeros(len(summary))
    for state in ["Inactive", "Active", "Intermediate"]:
        vals = summary[state].to_numpy()
        ax2.barh(y_pos, vals, left=left, color=palette[state], alpha=0.85, label=state)
        left += vals

    for idx, total_val in enumerate(left):
        ax2.text(total_val + 1.5, idx, f"{int(total_val)}", va="center", fontsize=8.5, fontweight="bold")

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(summary.index, fontsize=8.2)
    ax2.set_xlim(0, max(left) * 1.16)
    ax2.set_xlabel("Number of PDB Structures", fontsize=10)
    ax2.set_title("B. Construct Engineering vs. Conformational State", fontsize=10.5, fontweight="bold")
    ax2.legend(loc="lower right", fontsize=8)
    ax2.grid(True, axis="x", alpha=0.22)

    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
