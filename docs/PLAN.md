# Implementation Plan · ADORA2A Retrospective Audit

Reproducible case study on the human Adenosine $A_{2A}$ receptor (`ADORA2A`, UniProt `P29274`, ChEMBL `CHEMBL251`).  
**Objective:** demonstrate in a clean, fast, and verifiable repository the scientific rigor required to design, audit, and evaluate a retrospective GPCR virtual screening campaign—and communicate those findings clearly to an AI / data engineering team.

---

## Guiding Principles & Explicit *Go / No-Go* Criteria

- **No data without provenance:** every bioactivity record retains its `assay_chembl_id`, `standard_type`, `standard_relation`, and `confidence_score`.
- **Never mix apples and oranges:** strictly separate equilibrium binding assays ($K_i$) from $IC_{50}$ / functional cellular readouts, and distinguish agonists from antagonists.
- **Benchmark bias control:** a *dummy* classifier trained solely on 1D physicochemical properties (`MW`, `cLogP`, `TPSA`, `HBD`, `HBA`, `RotBonds`, `FormalCharge`) must yield a cross-validated `ROC-AUC` near `0.50` on the matched active + decoy benchmark. If it exceeds `~0.60`, the decoy set is biased and rejected.
- **Honest screening metrics:** always report **EF1%**, **EF5%**, **BEDROC ($\alpha = 20$)**, **ROC-AUC**, and **PR-AUC**—never *accuracy* as a primary metric.
- **Explicit failure reporting:** transparently document negative findings (e.g., lack of linear correlation between docking score and experimental $pK_i$, rigid-docking failure on dehydrated pockets, or performance degradation under Murcko scaffold splitting).

---

## Implementation Phases

### Phase 1 · GPCR Structural Audit (`structures.py`)
- [x] Query the **GPCRdb** REST API (`/services/structure/protein/aa2ar_human/`) and **RCSB PDB** GraphQL metadata.
- [x] Classify all 188 structures by conformational state (`Active`, `Intermediate`, `Inactive`), resolution, method (X-ray / cryo-EM), co-crystallized ligand modality (agonist vs. antagonist), fusion partners (`BRIL` / `T4L`), and thermostabilizing mutations (`StaR2`).
- [x] Generate comparative audit report in `reports/structure_audit.md` + summary figure (resolution vs. conformational state & construct class).
- [x] **Explicit decision:** select 1 reference structure for antagonists (inactive: `4EIY`) and 1 for agonists (active: `9EE8`), justifying every rejected structure.

### Phase 2 · Bioactivity Curation & Retrospective Set (`curation.py` + `decoys.py`)
- [x] Fetch `CHEMBL251` bioactivities from **ChEMBL** (with disk caching in `data/cache/` for offline reproducibility).
- [x] Filter by `confidence_score >= 8`, `standard_relation == '='`, normalized `nM` units, and strictly separate cohorts ($K_i$ vs. $IC_{50}$).
- [x] RDKit chemical standardization: largest organic fragment (salt stripping), neutralization, deduplication by `InChIKey` (computing median $pK_i$ and flagging $> 1.0$ log-unit inter-assay discrepancies).
- [x] Contextual structural alert auditing (PAINS/Brenk), verifying that clinically validated drugs and reference probes (`Caffeine`, `Istradefylline`, `ZM241385`, `Tozadenant`, `Vipadenant`, `Imaradenant`) are not blindly discarded.
- [x] Generate **7D property-matched decoys** (`MW`, `cLogP`, `TPSA`, `HBD`, `HBA`, `RotatableBonds`, `FormalCharge`) with low topological similarity (`ECFP4 Tanimoto < 0.35` to every active).
- [x] Verify with 1D physicochemical classifiers (`LogisticRegression` and `RandomForest`) that the benchmark passes the `< 0.60` ROC-AUC *go/no-go* gate.

### Phase 3 · Retrospective Benchmark & Honest Metrics (`metrics.py` + `benchmark.py`)
- [x] Implement and unit-test (`tests/test_metrics.py`) **EF1%**, **EF5%**, **BEDROC ($\alpha=20$)**, **ROC-AUC**, **PR-AUC**, and *accuracy* (contrasting why accuracy misleads at a 2–3% hit rate, and why global ROC-AUC hides early-recognition failures).
- [x] Run a reproducible 4-level comparison on the exact same 1,000-compound set:
  1. **1D Negative Control:** physicochemical classifier (`MW`, `cLogP`, `TPSA`, `Charge`, `HBD`, `HBA`, `RotBonds`) proving absence of trivial benchmark bias.
  2. **2D Baseline:** ECFP4 Tanimoto similarity against the co-crystallized antagonist `ZM241385` (plus non-analogue subset & **Murcko scaffold split** generalization audit).
  3. **Fast 3D Shape Baseline:** RDKit 3D conformers (`ETKDGv3` + `MMFF94`) scored via `USRCAT` and `CrippenO3A` Shape Tanimoto against crystallographic `ZM241385`.
  4. **Structure-Based Docking (AutoDock Vina):** prepare pH 7.4 wild-type `4EIY` with the orthosteric box centered on `ZM241385` and score all benchmark ligands; publish negative early-enrichment findings transparently alongside cognate redocking diagnostics.
- [x] Quantify **docking / similarity score vs. experimental $pK_i$** correlation across 60 measured binders spanning 5 orders of magnitude (*score $\neq$ affinity*).
- [x] Generate publication-quality figures in `reports/figures/`.

### Phase 4 · Executive Summary & Release (`reports/executive_summary.md` + `README.md`)
- [x] Write a 1-page executive summary tailored to AI/data engineers: therapeutic framing, *go/no-go* gates, curation choices, honest negative results, and production pipeline recommendations.
- [x] Complete `README.md` and executed walkthrough notebook (`notebooks/01_adora2a_retrospective_audit.ipynb`).
- [x] Verify `pytest` and `ruff check .` green + clean git commit.
