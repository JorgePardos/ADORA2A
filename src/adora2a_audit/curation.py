"""Bioactivity curation, RDKit chemical standardization, and contextual PAINS/Brenk auditing.

Enforces strict provenance preservation (`assay_chembl_id`, `standard_type`,
`standard_relation`, `confidence_score`), separates equilibrium binding (`Ki`)
from `IC50` / functional cohorts, deduplicates standardized parent structures by
`InChIKey`, and flags inter-assay discrepancies (> 1.0 log unit) without blind
deletion of reference ligands triggered by structural alert catalogs.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, FilterCatalog
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

_LARGEST_FRAG = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
_UNCHARGER = rdMolStandardize.Uncharger()


def _build_filter_catalog(catalog_enum: FilterCatalog.FilterCatalogParams.FilterCatalogs) -> FilterCatalog.FilterCatalog:
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(catalog_enum)
    return FilterCatalog.FilterCatalog(params)


_PAINS_CATALOG = _build_filter_catalog(
    FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS
)
_BRENK_CATALOG = _build_filter_catalog(
    FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK
)

# Canonical SMARTS hallmarks for known orthosteric A2A agonist chemotypes:
# 1) Ribofuranosyl nucleosides (adenosine, NECA, CGS21680, regadenoson, UK-432097)
# 2) Non-nucleoside 3,5-dicyanopyridine partial agonists (BAY 60-6583 / LUF5834)
_RIBOSE_AGONIST_SMARTS = Chem.MolFromSmarts("[#7,#6]1:[#6,#7]:[#6,#7]:[#6,#7]:[#6,#7]1[C@H,CH]2O[C@H,CH]([CH2,C])[C@@H,CH](O)[C@H,CH]2O")
_RIBOSE_GENERIC_SMARTS = Chem.MolFromSmarts("O1[C]([#7])[C](O)[C](O)[C]1")
_DICYANOPYRIDINE_AGONIST_SMARTS = Chem.MolFromSmarts("N#Cc1c(N)nc(S)c(C#N)c1")

# Well-characterized ADORA2A reference ligands to verify during curation & alert audit
REFERENCE_LIGANDS: dict[str, dict[str, str]] = {
    "ZM241385": {
        "smiles": "Oc1ccc(cc1)CCNc1nc(N)n2c(n1)nc(n2)c1ccco1",
        "role": "Antagonist / Inverse agonist (4EIY co-crystal probe)",
        "clinical_status": "Preclinical gold-standard probe",
    },
    "Istradefylline": {
        "smiles": "COc1ccc(/C=C/c2nc3c(c(=O)n(CC)c(=O)n3CC)n2C)cc1OC",
        "role": "Antagonist (FDA-approved for Parkinson's disease, Nourianz)",
        "clinical_status": "Approved (Phase 4)",
    },
    "Caffeine": {
        "smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        "role": "Antagonist (non-selective methylxanthine, 5MZP co-crystal)",
        "clinical_status": "Approved (Phase 4)",
    },
    "Theophylline": {
        "smiles": "Cn1c(=O)n(C)c2c(c1=O)[nH]cn2",
        "role": "Antagonist (methylxanthine, 5MZJ co-crystal)",
        "clinical_status": "Approved (Phase 4)",
    },
    "Preladenant": {
        "smiles": "COCCn1nc(-c2ccco2)nc2nc(NCCN3CCN(c4ccc5ncc(C)nc5c4)CC3)nc21",
        "role": "Antagonist (selective clinical candidate, SCH-420814)",
        "clinical_status": "Phase 3",
    },
    "Vipadenant": {
        "smiles": "Nc1nc(c2ccco2)c2c(n1)n(nn2)Cc1ccc(c(c1)C)N",
        "role": "Antagonist (5OLH co-crystal)",
        "clinical_status": "Phase 2",
    },
    "Tozadenant": {
        "smiles": "COc1ccc(c2c1nc(s2)NC(=O)N1CCC(CC1)(C)O)N1CCOCC1",
        "role": "Antagonist (5OLO co-crystal)",
        "clinical_status": "Phase 3",
    },
    "Imaradenant": {
        "smiles": "Nc1nnc(c(n1)c1ccc(cc1)F)c1cc(C)nc(c1)Cl",
        "role": "Antagonist (AZD4635, 6GT3 co-crystal, immuno-oncology)",
        "clinical_status": "Phase 2",
    },
    "Adenosine": {
        "smiles": "OC[C@H]1O[C@H]([C@@H]([C@@H]1O)O)n1cnc2c1ncnc2N",
        "role": "Endogenous full agonist (2YDO / 9EE8 co-crystal)",
        "clinical_status": "Approved (Phase 4)",
    },
    "NECA": {
        "smiles": "CCNC(=O)[C@H]1O[C@H]([C@@H]([C@@H]1O)O)n1cnc2c1ncnc2N",
        "role": "Agonist (2YDV / 5G53 co-crystal)",
        "clinical_status": "Preclinical probe",
    },
    "CGS21680": {
        "smiles": "CCNC(=O)[C@H]1O[C@H]([C@@H]([C@@H]1O)O)n1cnc2c1nc(NCCc1ccc(cc1)CCC(=O)O)nc2N",
        "role": "Selective A2A agonist (4UHR co-crystal)",
        "clinical_status": "Preclinical probe",
    },
    "Regadenoson": {
        "smiles": "CNC(=O)c1cnn(-c2nc(N)c3ncn([C@@H]4O[C@H](CO)[C@@H](O)[C@H]4O)c3n2)c1",
        "role": "Agonist (FDA-approved myocardial perfusion agent, Lexiscan)",
        "clinical_status": "Approved (Phase 4)",
    },
}


@dataclass(frozen=True)
class StandardizedCompound:
    """Result of RDKit salt-stripping, neutralization, and identifier generation."""

    canonical_smiles: str
    inchikey: str
    murcko_smiles: str
    mw: float
    pains_alerts: tuple[str, ...]
    brenk_alerts: tuple[str, ...]


def standardize_smiles(smiles: str) -> StandardizedCompound | None:
    """Strip salts, neutralize ionizable groups, and compute parent InChIKey + alerts."""
    if not smiles or not isinstance(smiles, str):
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        frag = _LARGEST_FRAG.choose(mol)
        uncharged = _UNCHARGER.uncharge(frag)
        Chem.SanitizeMol(uncharged)
    except (RuntimeError, ValueError):
        return None

    # Reject inorganic / tiny fragments (< 6 heavy atoms) or organometallics
    allowed_atomic_nums = {1, 5, 6, 7, 8, 9, 15, 16, 17, 35, 53}
    if uncharged.GetNumHeavyAtoms() < 6:
        return None
    if any(a.GetAtomicNum() not in allowed_atomic_nums for a in uncharged.GetAtoms()):
        return None

    can_smi = Chem.MolToSmiles(uncharged, isomericSmiles=True)
    inchikey = Chem.MolToInchiKey(uncharged)
    if not inchikey:
        return None

    try:
        scaf = MurckoScaffold.GetScaffoldForMol(uncharged)
        murcko_smi = Chem.MolToSmiles(scaf, isomericSmiles=False) if scaf.GetNumHeavyAtoms() > 0 else can_smi
    except (RuntimeError, ValueError):
        murcko_smi = can_smi

    pains = tuple(sorted({e.GetDescription() for e in _PAINS_CATALOG.GetMatches(uncharged)}))
    brenk = tuple(sorted({e.GetDescription() for e in _BRENK_CATALOG.GetMatches(uncharged)}))
    mw = float(Descriptors.MolWt(uncharged))

    return StandardizedCompound(
        canonical_smiles=can_smi,
        inchikey=inchikey,
        murcko_smiles=murcko_smi,
        mw=mw,
        pains_alerts=pains,
        brenk_alerts=brenk,
    )


def build_reference_inchikey_map() -> dict[str, tuple[str, str, str]]:
    """Map standardized InChIKey -> (name, role, clinical_status) for reference ligands."""
    mapping: dict[str, tuple[str, str, str]] = {}
    for name, info in REFERENCE_LIGANDS.items():
        std = standardize_smiles(info["smiles"])
        if std is not None:
            # Match by 14-char connectivity layer as well in case stereoisomer annotation varies
            mapping[std.inchikey] = (name, info["role"], info["clinical_status"])
            mapping[std.inchikey[:14]] = (name, info["role"], info["clinical_status"])
    return mapping


def _build_known_modality_lookup(cache_dir: Path) -> dict[str, str]:
    """Build connectivity-InChIKey -> modality ('Agonist' | 'Antagonist') from ChEMBL + GPCRdb."""
    lookup: dict[str, str] = {}

    gpcr_path = cache_dir / "gpcrdb_aa2ar_human.json"
    if gpcr_path.exists():
        gpcr = json.loads(gpcr_path.read_text(encoding="utf-8"))
        for entry in gpcr:
            for lig in entry.get("ligands", []):
                smi = lig.get("SMILES")
                func = (lig.get("function") or "").lower()
                if not smi or not func:
                    continue
                std = standardize_smiles(smi)
                if std is None:
                    continue
                if "agonist" in func and "inverse" not in func:
                    lookup[std.inchikey[:14]] = "Agonist"
                elif "antagonist" in func or "inverse" in func:
                    lookup[std.inchikey[:14]] = "Antagonist"

    for info in REFERENCE_LIGANDS.values():
        std = standardize_smiles(info["smiles"])
        if std is not None:
            role_low = info["role"].lower()
            lookup[std.inchikey[:14]] = (
                "Agonist" if ("agonist" in role_low and "inverse" not in role_low and "antagonist" not in role_low) else "Antagonist"
            )
    return lookup


def infer_modality(
    canonical_smiles: str,
    inchikey: str,
    molecule_chembl_ids: list[str],
    assay_descriptions: list[str],
    mech_map: dict[str, str],
    struct_modality_map: dict[str, str],
) -> str:
    """Classify ligand functional modality ('Antagonist', 'Agonist', or 'Binding_Unclassified').

    Priority:
    1. Curated ChEMBL mechanism of action or GPCRdb co-crystal annotation.
    2. Explicit functional / pharmacological annotation in ChEMBL assay description.
    3. Canonical A2A agonist chemotype check (ribofuranosyl nucleoside or 3,5-dicyanopyridine).
    4. Radioligand antagonist displacement context -> 'Antagonist (radioligand displacement)'.
    """
    for cid in molecule_chembl_ids:
        if cid in mech_map:
            act = mech_map[cid].upper()
            if "AGONIST" in act and "INVERSE" not in act and "ANTAGONIST" not in act:
                return "Agonist"
            if "ANTAGONIST" in act or "INVERSE" in act:
                return "Antagonist"

    conn = inchikey[:14]
    if conn in struct_modality_map:
        return struct_modality_map[conn]

    mol = Chem.MolFromSmiles(canonical_smiles)
    if mol is not None and (
        (_RIBOSE_AGONIST_SMARTS and mol.HasSubstructMatch(_RIBOSE_AGONIST_SMARTS))
        or (_RIBOSE_GENERIC_SMARTS and mol.HasSubstructMatch(_RIBOSE_GENERIC_SMARTS))
        or (_DICYANOPYRIDINE_AGONIST_SMARTS and mol.HasSubstructMatch(_DICYANOPYRIDINE_AGONIST_SMARTS))
    ):
        return "Agonist"

    joined_desc = " || ".join(assay_descriptions).lower()
    if re.search(r"\b(antagonist\s+activity|inverse\s+agonist|antagonism)\b", joined_desc):
        return "Antagonist"
    if re.search(r"\b(agonist\s+activity|full\s+agonist|partial\s+agonist)\b", joined_desc):
        return "Agonist"
    if re.search(r"zm[- ]?241385|sch[- ]?58261|istradefylline|xac|msx[- ]?2", joined_desc):
        return "Antagonist (binding)"
    return "Binding_Unclassified"


def curate_chembl_activities(
    cache_dir: str | Path,
    min_confidence: int = 8,
    active_pki_threshold: float = 7.0,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Load, filter, standardize, and deduplicate CHEMBL251 bioactivity records.

    Returns:
        raw_filtered_df: Measurement-level records with full provenance (`assay_chembl_id`, etc.)
        dedup_ki_df: Compound-level deduplicated equilibrium `Ki` binding dataset (`assay_type == 'B'`)
        waterfall: Step-by-step audit counts documenting every excluded record
    """
    cache_path = Path(cache_dir)
    acts = json.loads((cache_path / "chembl251_activities_ki_ic50.json").read_text(encoding="utf-8"))
    assays = {
        a["assay_chembl_id"]: a
        for a in json.loads((cache_path / "chembl251_assays.json").read_text(encoding="utf-8"))
    }
    mechs_raw = json.loads((cache_path / "chembl251_mechanisms.json").read_text(encoding="utf-8"))
    mech_map: dict[str, str] = {}
    for m in mechs_raw:
        action = m.get("action_type") or ""
        for k in ("molecule_chembl_id", "parent_molecule_chembl_id"):
            if m.get(k):
                mech_map[m[k]] = action

    struct_modality_map = _build_known_modality_lookup(cache_path)
    ref_map = build_reference_inchikey_map()

    df = pd.DataFrame(acts)
    waterfall: dict[str, int] = {"01_raw_ki_ic50_records": len(df)}

    # Attach assay provenance
    df["confidence_score"] = df["assay_chembl_id"].map(
        lambda x: assays.get(x, {}).get("confidence_score")
    )
    df["assay_type"] = df["assay_chembl_id"].map(lambda x: assays.get(x, {}).get("assay_type"))
    df["assay_description"] = df["assay_chembl_id"].map(
        lambda x: assays.get(x, {}).get("description") or ""
    )

    # Filter 1: confidence_score >= min_confidence
    df = df[df["confidence_score"].fillna(0).astype(int) >= min_confidence].copy()
    waterfall["02_confidence_ge_8"] = len(df)

    # Filter 2: exact relation '=' (exclude censored '>' / '<' values from quantitative medians)
    df = df[df["standard_relation"] == "="].copy()
    waterfall["03_exact_relation_eq"] = len(df)

    # Filter 3: valid positive nanomolar concentration & SMILES
    df["value_nm"] = pd.to_numeric(df["normalized_value_nM"], errors="coerce")
    df = df[(df["value_nm"] > 0) & df["canonical_smiles"].notna()].copy()
    df["p_activity"] = 9.0 - df["value_nm"].apply(math.log10)
    waterfall["04_valid_nm_and_smiles"] = len(df)

    # Filter 4: drop explicit data_validity_comment flags from ChEMBL ('Outside typical range', etc.)
    df = df[df["data_validity_comment"].isna()].copy()
    waterfall["05_no_validity_warning"] = len(df)

    # Chemical standardization cache (avoid recomputing identical SMILES)
    unique_smiles = df["canonical_smiles"].unique().tolist()
    std_cache = {smi: standardize_smiles(smi) for smi in unique_smiles}
    df["std_obj"] = df["canonical_smiles"].map(std_cache)
    df = df[df["std_obj"].notna()].copy()

    df["std_smiles"] = df["std_obj"].map(lambda s: s.canonical_smiles)
    df["inchikey"] = df["std_obj"].map(lambda s: s.inchikey)
    df["murcko_smiles"] = df["std_obj"].map(lambda s: s.murcko_smiles)
    df["mw"] = df["std_obj"].map(lambda s: s.mw)
    df["pains_alerts"] = df["std_obj"].map(lambda s: "; ".join(s.pains_alerts))
    df["brenk_alerts"] = df["std_obj"].map(lambda s: "; ".join(s.brenk_alerts))
    df = df.drop(columns=["std_obj"])
    waterfall["06_rdkit_standardized_organic"] = len(df)

    # Strictly separate cohorts: Equilibrium Binding Ki (B) vs Binding IC50 (B) vs Functional IC50 (F)
    df["cohort"] = (
        df["standard_type"].astype(str) + "_" + df["assay_type"].fillna("U").astype(str)
    )
    waterfall["07_cohort_Ki_B_measurements"] = int((df["cohort"] == "Ki_B").sum())
    waterfall["07b_cohort_IC50_B_measurements"] = int((df["cohort"] == "IC50_B").sum())
    waterfall["07c_cohort_IC50_F_measurements"] = int((df["cohort"] == "IC50_F").sum())

    # Deduplicate primary Ki_B cohort by parent InChIKey while preserving ALL provenance IDs
    ki_df = df[df["cohort"] == "Ki_B"].copy()
    dedup_rows = []
    for ik, grp in ki_df.groupby("inchikey", sort=True):
        pvals = grp["p_activity"].to_numpy(dtype=float)
        pki_med = float(pd.Series(pvals).median())
        pki_iqr = float(pd.Series(pvals).quantile(0.75) - pd.Series(pvals).quantile(0.25))
        pki_span = float(pvals.max() - pvals.min())
        mol_ids = sorted(grp["molecule_chembl_id"].dropna().unique().tolist())
        assay_ids = sorted(grp["assay_chembl_id"].dropna().unique().tolist())
        doc_ids = sorted(grp["document_chembl_id"].dropna().unique().tolist())
        descs = grp["assay_description"].dropna().unique().tolist()
        std_smi = grp["std_smiles"].iloc[0]

        modality = infer_modality(
            canonical_smiles=std_smi,
            inchikey=ik,
            molecule_chembl_ids=mol_ids,
            assay_descriptions=descs,
            mech_map=mech_map,
            struct_modality_map=struct_modality_map,
        )
        ref_hit = ref_map.get(ik) or ref_map.get(ik[:14])
        pains_str = grp["pains_alerts"].iloc[0]
        brenk_str = grp["brenk_alerts"].iloc[0]

        dedup_rows.append(
            {
                "inchikey": ik,
                "canonical_smiles": std_smi,
                "murcko_smiles": grp["murcko_smiles"].iloc[0],
                "primary_chembl_id": mol_ids[0] if mol_ids else "",
                "molecule_chembl_ids": ";".join(mol_ids),
                "assay_chembl_ids": ";".join(assay_ids),
                "document_chembl_ids": ";".join(doc_ids),
                "standard_type": "Ki",
                "standard_relation": "=",
                "min_confidence_score": int(grp["confidence_score"].min()),
                "max_confidence_score": int(grp["confidence_score"].max()),
                "n_measurements": len(grp),
                "n_distinct_assays": len(assay_ids),
                "pki_median": round(pki_med, 4),
                "ki_median_nm": round(10 ** (9.0 - pki_med), 3),
                "pki_iqr": round(pki_iqr, 4),
                "pki_span": round(pki_span, 4),
                "inter_assay_discrepancy": bool(pki_span > 1.0),
                "modality": modality,
                "mw": round(float(grp["mw"].iloc[0]), 2),
                "has_pains": bool(pains_str),
                "pains_alerts": pains_str,
                "has_brenk": bool(brenk_str),
                "brenk_alerts": brenk_str,
                "reference_name": ref_hit[0] if ref_hit else "",
                "reference_role": ref_hit[1] if ref_hit else "",
                "is_potent_active": bool(
                    pki_med >= active_pki_threshold and modality != "Agonist" and pki_span <= 1.5
                ),
            }
        )

    dedup_ki_df = pd.DataFrame(dedup_rows).sort_values("pki_median", ascending=False).reset_index(drop=True)
    waterfall["08_unique_inchikey_Ki_B"] = len(dedup_ki_df)
    waterfall["09_discrepant_gt_1log_Ki_B"] = int(dedup_ki_df["inter_assay_discrepancy"].sum())
    waterfall["10_agonists_excluded_from_inactive_screen"] = int(
        (dedup_ki_df["modality"] == "Agonist").sum()
    )
    waterfall["11_potent_non_agonist_actives_pKi_ge_7"] = int(
        dedup_ki_df["is_potent_active"].sum()
    )

    provenance_cols = [
        "activity_id",
        "molecule_chembl_id",
        "assay_chembl_id",
        "document_chembl_id",
        "standard_type",
        "standard_relation",
        "confidence_score",
        "assay_type",
        "cohort",
        "value_nm",
        "p_activity",
        "inchikey",
        "std_smiles",
        "murcko_smiles",
        "pains_alerts",
        "brenk_alerts",
    ]
    existing_cols = [c for c in provenance_cols if c in df.columns]
    return df[existing_cols].reset_index(drop=True), dedup_ki_df, waterfall


def audit_reference_ligand_alerts() -> pd.DataFrame:
    """Evaluate PAINS and Brenk structural alerts on canonical ADORA2A reference ligands."""
    rows = []
    for name, info in REFERENCE_LIGANDS.items():
        std = standardize_smiles(info["smiles"])
        if std is None:
            continue
        rows.append(
            {
                "name": name,
                "role": info["role"],
                "clinical_status": info["clinical_status"],
                "inchikey": std.inchikey,
                "mw": round(std.mw, 1),
                "has_pains": len(std.pains_alerts) > 0,
                "pains_alerts": "; ".join(std.pains_alerts) or "None",
                "has_brenk": len(std.brenk_alerts) > 0,
                "brenk_alerts": "; ".join(std.brenk_alerts) or "None",
            }
        )
    return pd.DataFrame(rows)
