# ADORA2A Retrospective Virtual Screening & Curation Audit

A reproducible case study on the human **Adenosine $A_{2A}$ receptor** (`ADORA2A`, UniProt [`P29274`](https://www.uniprot.org/uniprotkb/P29274), ChEMBL [`CHEMBL251`](https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL251/)) demonstrating GPCR structure selection, provenance-preserving bioactivity curation, physicochemical decoy-bias control, and honest retrospective evaluation across 1D, 2D, 3D ligand-shape, and structure-based docking (**AutoDock Vina v1.2.7**) pipelines.

---

## Why This Project Exists

Virtual screening campaigns and machine-learning potency models frequently fail prospectively because of four avoidable upstream artefacts:

1. **Blind structure selection (`ORDER BY resolution ASC`):** docking into incompatible conformational states (agonist vs. antagonist) or engineered constructs carrying thermostabilizing mutations directly inside the orthosteric pocket (`StaR2` `T88A`/`S277A`).
2. **Uncurated bioactivity data without provenance:** mixing incompatible assay readouts (equilibrium binding $K_i$ vs. functional cellular $IC_{50}$), dropping `assay_chembl_id`, or applying PAINS/Brenk filters blindly and discarding clinically validated drugs.
3. **Trivial physicochemical benchmark leakage:** evaluating against unmatched decoys where a trivial 1D model (`MW`, `cLogP`, `TPSA`) achieves $\text{ROC-AUC} > 0.80$ without inspecting 2D chemical topology or 3D pocket geometry.
4. **Misleading metrics & confusing screening score with binding affinity:** reporting classification *accuracy* (~97.5%) or global `ROC-AUC` while ignoring early enrichment (**EF1%**, **BEDROC**), or assuming linear correlation between docking scores and experimental $pK_i$.

---

## Key Empirical Findings & Figures

### 1. Audit of All 188 ADORA2A Structures in PDB / GPCRdb ([reports/structure_audit.md](reports/structure_audit.md))

- **66.0% (`124 / 188`)** of human $A_{2A}$ PDB structures rely on the thermostabilized `StaR2` construct—including the two highest-resolution entries (`5NM4` at 1.70 Å and `5IU4` at 1.72 Å)—which mutate **`Thr88Ala` ($3.36$)** and **`Ser277Ala` ($7.42$)** directly inside the orthosteric binding pocket.
- **0 of 15 active-state structures (`0%`)** retain the conserved Class A allosteric $\text{Na}^+$ ion at $\text{Asp52}^{2.50}$, compared with **84.3% (`145 / 172`)** of inactive structures (**100%** at $\le 1.90\text{ Å}$).
- **Explicit reference choice:** **`4EIY` (1.80 Å, 100% wild-type 7TM pocket + allosteric $\text{Na}^+$ + `ZM241385`)** is selected for antagonist screening, and **`9EE8` (2.63 Å, $G_s$-coupled cryo-EM complex with endogenous adenosine)** for agonist studies.

![ADORA2A Structural Audit](reports/figures/structure_audit_resolution_state.png)

---

### 2. `CHEMBL251` Curation & 1D Decoy-Bias *Go / No-Go* Gate

- Starting from **12,105** raw `Ki`/`IC50` records in ChEMBL, exact measurements (`confidence_score >= 8`, `standard_relation == '='`, normalized `nM`, no validity warnings) are separated by assay type, yielding **`5,901` equilibrium binding `Ki_B` measurements $\rightarrow$ `4,972` unique parent structures (`InChIKey`)**, flagging `44` compounds with $> 1$ log-unit inter-assay discrepancy and excluding `546` agonists from the inactive-state screen.
- **Contextual structural alerts:** blind PAINS/Brenk filtering would have silently discarded three co-crystallized Phase 2/3 clinical $A_{2A}$ antagonists (`Tozadenant`, `Vipadenant`, and `Imaradenant`).
- **1D leakage gate passed (*GO*):** against random ChEMBL decoys, a 1D classifier (`MW`, `cLogP`, `TPSA`, `HBD`, `HBA`, `RotBonds`, `Charge`) reaches **ROC-AUC = `0.81` (LogReg) / `0.78` (RF)** driven by the elevated `TPSA` of nitrogen-rich $A_{2A}$ heterocycles. After multivariate 7D matching + topological dissimilarity filtering (`ECFP4 < 0.35`), 1D ROC-AUC drops to **`0.48` / `0.35`** ($< 0.60$).

![1D Decoy Bias Audit](reports/figures/decoy_property_matching_audit.png)

---

### 3. Four-Level Retrospective Benchmark (1,000 Compounds · 2.5% Hit Rate)

| Screening Method | EF1% (95% CI) | EF5% | BEDROC ($\alpha=20$) | ROC-AUC | PR-AUC | Accuracy @ Top 1% |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **1. 1D Physicochemical (Control)** | `0.0x` `(0.0–0.0)` | `0.0x` | `0.026` | `0.346` | `0.019` | `96.5%` |
| **2. 2D ECFP4 (vs. `ZM241385`)** | **`36.0x`** `(16.0–40.0)` | **`9.6x`** | **`0.519`** | `0.725` | **`0.433`** | `98.3%` |
| **3. 3D Shape + USRCAT (`ETKDGv3` + `O3A`)** | `8.0x` `(0.0–20.0)` | `6.4x` | `0.355` | `0.722` | `0.166` | `96.9%` |
| **4. 3D Docking (`AutoDock Vina v1.2.7` · `4EIY`)** | **`0.0x`** `(0.0–0.0)` | `0.8x` | `0.063` | `0.584` | `0.033` | `96.5%` |

- **Why `ROC-AUC` and `Accuracy` mislead:** `2D ECFP4` (`0.725`) and `3D Shape + USRCAT` (`0.722`) have virtually identical global `ROC-AUC`, yet `2D ECFP4` recovers **4.5× more true actives in the Top 1% (`EF1% = 36.0x` vs. `8.0x`)**. Meanwhile, predicting every molecule inactive yields **`97.5%` accuracy**.
- **Transparent negative result on dehydrated rigid-receptor docking:** once 1D molecular-weight and `TPSA` shortcuts are removed, rigid water-stripped `AutoDock Vina` drops to `EF1% = 0.0x` (`ROC-AUC = 0.584`; rising to `0.625` on the non-`ZM241385`-analogue subset). Cognate redocking of `ZM241385` pinpoints the physical cause: removing bridging waters `HOH 2503–2505` makes an inverted pose (`-9.08 kcal/mol`) nearly isoenergetic (`ΔE = 0.15 kcal/mol`) with the native crystal orientation (`-8.93 kcal/mol`, furan oxygen within `0.41 Å`).

![ROC and Early Enrichment Curves](reports/figures/benchmark_roc_and_enrichment.png)

![Accuracy and ROC-AUC vs Early Enrichment Metrics](reports/figures/metrics_accuracy_vs_enrichment.png)

---

### 4. Screening Score $\neq$ Experimental Binding Affinity ($pK_i$)

Across **60 experimentally measured $K_i$ binders** spanning more than 5 orders of magnitude ($pK_i \in [4.0,\; 11.0]$), neither AutoDock Vina binding energy ($r = +0.16$, $R^2 = 0.03$, $p = 0.22$) nor 2D ECFP4 similarity ($r = +0.25$, $R^2 = 0.06$) correlates linearly with experimental $pK_i$.

![Screening Score vs Experimental pKi](reports/figures/affinity_vs_score_scatter.png)

---

## Repository Structure

- [`src/adora2a_audit/`](src/adora2a_audit/)
  - [`structures.py`](src/adora2a_audit/structures.py) — combined GPCRdb + RCSB PDB audit and atomic pocket verification with `gemmi`.
  - [`curation.py`](src/adora2a_audit/curation.py) — provenance-preserving `CHEMBL251` curation, RDKit standardization, and contextual PAINS/Brenk auditing.
  - [`decoys.py`](src/adora2a_audit/decoys.py) — 7D property-matched decoy selection and 1D physicochemical leakage *go/no-go* gate.
  - [`metrics.py`](src/adora2a_audit/metrics.py) — tie-aware analytical implementations of `EF1%`, `EF5%`, `BEDROC` ($\alpha=20$), `ROC-AUC`, `PR-AUC`, and stratified bootstrap CIs.
  - [`benchmark.py`](src/adora2a_audit/benchmark.py) — pH 7.4 `4EIY` receptor preparation, cognate redocking audit, 2D similarity + *Murcko scaffold split*, 3D shape (`USRCAT` + `CrippenO3A`), and parallel `--batch` AutoDock Vina docking.
- [`scripts/run_audit_pipeline.py`](scripts/run_audit_pipeline.py) — deterministic end-to-end CLI pipeline.
- [`notebooks/01_adora2a_retrospective_audit.ipynb`](notebooks/01_adora2a_retrospective_audit.ipynb) — executed walkthrough notebook with rendered tables and figures.
- [`reports/structure_audit.md`](reports/structure_audit.md) — detailed audit across all 188 experimental `ADORA2A` structures.
- [`reports/executive_summary.md`](reports/executive_summary.md) — 1-page executive summary & production architecture recommendations for AI/data engineering teams.
- [`tests/`](tests/) — `pytest` unit test suite covering curation, metrics, structural checks, and decoy dissimilarity gates.

---

## Quick Reproduction

```bash
# 1. Install package in editable mode with development dependencies
pip install -e ".[dev]"

# 2. Run the end-to-end audit pipeline (reuses cached data in data/cache/ and processed scores)
python scripts/run_audit_pipeline.py

# 3. Run unit tests and linter
pytest -v
ruff check .
```

## Data Sources & Licensing

- **GPCRdb:** receptor state and co-crystallized ligand annotations retrieved via the [GPCRdb REST API](https://gpcrdb.org/).
- **RCSB PDB:** crystallographic metadata and atomic coordinates subject to the [RCSB PDB Usage Policy](https://www.rcsb.org/pages/usage-policy).
- **EMBL-EBI ChEMBL:** bioactivity records and chemical structures licensed under [CC BY-SA 3.0](https://chembl.gitbook.io/chembl-interface-documentation/about).
