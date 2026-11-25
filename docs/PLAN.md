# Implementation Plan · ADORA2A Retrospective Audit (v2)

Reproducible case study on the human Adenosine $A_{2A}$ receptor (`ADORA2A`, UniProt `P29274`, ChEMBL `CHEMBL251`).  
**Objective:** provide a transparent, verifiable audit of GPCR structure selection, bioactivity curation, decoy-bias control, and retrospective virtual screening—including a single-variable A/B test of crystallographic pocket hydration (`4EIY_dry` vs. `4EIY_hydrated`)—tailored for AI and data engineering teams.

---

## Guiding Principles & Explicit *Go / No-Go* Criteria

- **No data without provenance:** every bioactivity record retains its `assay_chembl_id`, `standard_type`, `standard_relation`, and `confidence_score`.
- **Strict assay & modality separation:** separate equilibrium binding assays ($K_i$) from $IC_{50}$ / functional cellular readouts, and exclude agonists from inactive-state antagonist screening.
- **Benchmark bias control:** 1D physicochemical classifiers (`MW`, `cLogP`, `TPSA`, `HBD`, `HBA`, `RotBonds`, `FormalCharge`) must yield cross-validated `ROC-AUC` $\le 0.60$ on the matched active + decoy benchmark.
- **Early-enrichment metrics:** report **EF1%**, **EF5%**, **BEDROC ($\alpha = 20$)**, **ROC-AUC**, and **PR-AUC** with stratified bootstrap confidence intervals—never *accuracy* as a primary metric.
- **Empirical discipline:** distinguish tested mechanisms from untested hypotheses, validate cognate redocking across multiple random seeds rather than a single seed, and report negative or modest A/B outcomes as measured.

---

## Implementation Phases

### Phase 1 · GPCR Structural & Pocket-Water Audit (`structures.py`)
- [x] Query the **GPCRdb** REST API (`/services/structure/protein/aa2ar_human/`) and **RCSB PDB** GraphQL metadata across all 188 entries.
- [x] Classify structures by conformational state (`Active`, `Intermediate`, `Inactive`), resolution, experimental method, ligand modality, intracellular fusion partners (`BRIL` / `T4L`), and `StaR2` thermostabilizing mutations (`T88A`, `S277A`).
- [x] Select `4EIY` (1.80 Å, wild-type 7TM pocket + resolved allosteric $\text{Na}^+$) as the inactive antagonist reference and `9EE8` (2.63 Å, $G_s$-coupled cryo-EM with adenosine) as the active reference.
- [x] Inspect crystallographic waters in `4EIY` with `gemmi` (`inspect_pocket_waters`) and document an objective four-criterion rule (`d_ZMA <= 4.5 Å`, `B < 30.0 Å²`, `occupancy == 1.0`, $\ge 1$ H-bond $\le 3.5\text{ Å}$) retaining 10 ordered pocket waters for the `4EIY_hydrated` arm.

### Phase 2 · Bioactivity Curation & 7D-Matched Decoy Set (`curation.py` + `decoys.py`)
- [x] Fetch and cache `CHEMBL251` bioactivities, assays, and mechanisms of action from **ChEMBL**.
- [x] Filter by `confidence_score >= 8`, `standard_relation == '='`, normalized `nM` units, and separate equilibrium binding `Ki_B` from `IC50` cohorts.
- [x] Standardize chemical structures with RDKit (salt stripping, neutralization, `InChIKey` deduplication, median $pK_i$, flagging $> 1.0$ log-unit inter-assay discrepancies).
- [x] Audit PAINS/Brenk alerts in pharmacological context (`Caffeine`, `Istradefylline`, `ZM241385`, `Tozadenant`, `Vipadenant`, `Imaradenant`).
- [x] Generate **7D property-matched decoys** (`MW`, `cLogP`, `TPSA`, `HBD`, `HBA`, `RotatableBonds`, `FormalCharge`) with `ECFP4 Tanimoto < 0.35` to all actives and verify `< 0.60` 1D `ROC-AUC`.

### Phase 3 · 10-Seed Cognate Redocking & 5-Series Retrospective Benchmark (`metrics.py` + `benchmark.py`)
- [x] Implement and unit-test (`tests/test_metrics.py`) **EF1%**, **EF5%**, **BEDROC ($\alpha=20$)**, **ROC-AUC**, **PR-AUC**, and *accuracy*.
- [x] Run **10 independent random seeds** of cognate `ZM241385` redocking on both `4EIY_dry` and `4EIY_hydrated`, reporting Pose-1 and best-pose symmetry-corrected heavy-atom RMSD, `< 2.0 Å` recovery rate, and $\Delta E$ (`native - non-native`).
- [x] Evaluate on the exact same 1,000-compound benchmark (`25` actives + `975` 7D-matched decoys, `2.5%` hit rate):
  1. **1D Physicochemical Control** (`RandomForest` & `LogisticRegression` OOF)
  2. **2D ECFP4 Similarity** (vs. `ZM241385`, plus non-analogue subset & **Murcko scaffold split** audit)
  3. **3D Ligand Shape + USRCAT** (`ETKDGv3` + `MMFF94` + `CrippenO3A` vs. crystallographic `ZMA`)
  4. **3D Structure-Based Docking — Dry Arm** (`AutoDock Vina v1.2.7` on `4EIY_dry`)
  5. **3D Structure-Based Docking — Hydrated Arm** (`AutoDock Vina v1.2.7` on `4EIY_hydrated`)
- [x] Evaluate score vs. experimental $pK_i$ correlation across 60 measured binders (`pKi` 4.0–11.0) for both dry and hydrated receptor arms.
- [x] Generate updated figures in `reports/figures/`.

### Phase 4 · Sober Documentation, Executed Notebook & Release (`reports/executive_summary.md` + `README.md`)
- [x] Update `reports/executive_summary.md` and `README.md` with data-driven mechanistic conclusions and explicit methodological notes.
- [x] Update and execute `notebooks/01_adora2a_retrospective_audit.ipynb`.
- [x] Verify `pytest -v` and `ruff check .` green + clean git commit on `main`.
