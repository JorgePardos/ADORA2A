# ADORA2A Retrospective Virtual Screening & Curation Audit

A reproducible case study on the human **Adenosine $A_{2A}$ receptor** (`ADORA2A`, UniProt [`P29274`](https://www.uniprot.org/uniprotkb/P29274), ChEMBL [`CHEMBL251`](https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL251/)) examining GPCR structure selection, provenance-preserving bioactivity curation, physicochemical decoy-bias control, and retrospective screening evaluation across 1D, 2D, 3D ligand-shape, and structure-based docking (**AutoDock Vina v1.2.7**, comparing dry vs. crystallographically hydrated `4EIY` receptor arms).

---

## Motivation

Virtual screening campaigns and potency models often fail prospectively due to four upstream issues:

1. **Unfiltered structure selection (`ORDER BY resolution ASC`):** docking into an incompatible conformational state (agonist vs. antagonist) or into thermostabilized constructs carrying point mutations inside the orthosteric pocket (`StaR2` `T88A`/`S277A`).
2. **Uncurated bioactivity records:** pooling incompatible assay readouts (equilibrium binding $K_i$ vs. functional cellular $IC_{50}$), dropping `assay_chembl_id`, or applying PAINS/Brenk catalogs as blind hard filters that remove validated clinical candidates.
3. **Physicochemical benchmark leakage:** evaluating against unmatched decoys where a 1D descriptor baseline (`MW`, `cLogP`, `TPSA`) achieves $\text{ROC-AUC} > 0.78$ without using 2D topology or 3D pocket geometry.
4. **Inadequate evaluation metrics & over-interpreting docking scores:** reporting classification *accuracy* (~97.5% at a 2.5% hit rate) or global `ROC-AUC` without early enrichment (**EF1%**, **BEDROC**), or treating empirical docking scores as proxies for experimental $pK_i$.

---

## Summary of Empirical Findings

### 1. Structural & Pocket-Water Audit Across 188 PDB Entries ([reports/structure_audit.md](reports/structure_audit.md))

- **66.0% (`124 / 188`)** of human $A_{2A}$ structures in the PDB use the thermostabilized `StaR2` construct—including the two highest-resolution entries (`5NM4` at 1.70 Å and `5IU4` at 1.72 Å)—which introduce **`Thr88Ala` ($3.36$)** and **`Ser277Ala` ($7.42$)** inside the orthosteric pocket.
- **0 of 15 active-state structures (`0%`)** retain the allosteric $\text{Na}^+$ ion at $\text{Asp52}^{2.50}$, compared with **84.3% (`145 / 172`)** of inactive structures (**100%** at $\le 1.90\text{ Å}$).
- **Reference structures:** **`4EIY` (1.80 Å, wild-type 7TM pocket + allosteric $\text{Na}^+$ + `ZM241385`)** is selected for antagonist screening, and **`9EE8` (2.63 Å, $G_s$-coupled cryo-EM complex with adenosine)** as the active-state reference.
- **Receptor preparation & A/B hydration protocol:** receptor PDBQT files (`4EIY_dry` and `4EIY_hydrated`) were prepared at pH 7.4 via a custom `RDKit`/`gemmi` writer with Gasteiger partial charges (noting that AutoDock Vina's `vina` scoring function ignores partial charges and uses only AutoDock atom types plus covalent `HD` connectivity). Inspecting `4EIY` with `gemmi` identifies 15 waters within `5.0 Å` of `ZM241385`; applying an objective four-part rule ($d_{\text{ZMA}} \le 4.5\text{ Å}$, $B < 30.0\text{ Å}^2$, occupancy $= 1.00$, $\ge 1$ H-bond $\le 3.5\text{ Å}$) retains **10 ordered waters** (`HOH 2517, 2518, 2521, 2524, 2525, 2527, 2531, 2541, 2583, 2585`) as rigid `OA`+`HD` (`XS_TYPE_O_DA`) receptor atoms in `4EIY_hydrated`.

![ADORA2A Structural Audit](reports/figures/structure_audit_resolution_state.png)

---

### 2. `CHEMBL251` Curation & 1D Decoy-Bias Gate

- From **12,105** raw `Ki`/`IC50` records in ChEMBL, filtering exact measurements (`confidence_score >= 8`, `standard_relation == '='`, normalized `nM`, no validity warnings) and isolating equilibrium binding assays yields **`5,901` `Ki_B` measurements $\rightarrow$ `4,972` unique parent compounds (`InChIKey`)**, flagging `44` compounds with $> 1$ log-unit inter-assay spread and excluding `546` agonists from the inactive-state screen.
- **Contextual structural alerts:** blind PAINS/Brenk deletion would remove three co-crystallized Phase 2/3 clinical $A_{2A}$ antagonists (`Tozadenant`, `Vipadenant`, and `Imaradenant`).
- **1D leakage audit (*GO*):** against random ChEMBL decoys, a 5-fold CV classifier using only 1D properties (`MW`, `cLogP`, `TPSA`, `HBD`, `HBA`, `RotBonds`, `Charge`) reaches **ROC-AUC = `0.81` (LogReg) / `0.78` (RF)** due to higher `TPSA` in $A_{2A}$ heterocycles. After 7D property matching + topological filtering (`ECFP4 < 0.35`), 1D ROC-AUC drops to **`0.484` (LogReg)** and **`0.346` (RF)**—passing the `< 0.60` gate, while noting that `0.346` in RF reflects slight inverted out-of-fold signal and small-sample variance across 25 locally matched actives rather than exact `0.50` uniform noise.

![1D Decoy Bias Audit](reports/figures/decoy_property_matching_audit.png)

---

### 3. 10-Seed Cognate Redocking & 1,000-Compound A/B Benchmark (2.5% Hit Rate)

#### 3.1 Multi-Seed Cognate Redocking (`ZM241385`, 10 Independent Random Seeds)

| Receptor Arm | Pose 1 RMSD (Å) | Pose 1 $< 2.0\text{ Å}$ | Best-of-9 RMSD (Å) | Any Pose $< 2.0\text{ Å}$ | $\Delta E = E_{\text{native}} - E_{\text{non-native}}$ (kcal/mol) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`4EIY_dry`** | `9.68 ± 0.02` | `0 / 10` | `1.34 ± 0.49` | `9 / 10` | `+0.188 ± 0.084` *(inverted pose favored)* |
| **`4EIY_hydrated` (10 waters)** | **`1.05 ± 0.15`** | **`10 / 10`** | **`1.05 ± 0.15`** | **`10 / 10`** | **`-0.254 ± 0.044`** *(native pose favored)* |

Across 10 seeds, `4EIY_dry` samples the native crystallographic geometry in 9/10 runs (`1.34 Å`) but consistently ranks a head-to-tail inverted pose first (`+0.19 kcal/mol` scoring inversion). Retaining the 10 ordered first-shell waters in `4EIY_hydrated` restores top-1 native recovery in **10/10 seeds (`1.05 ± 0.15 Å`)**.

#### 3.2 Retrospective Screening Performance (25 Actives + 975 7D-Matched Decoys)

| Screening Method | EF1% (95% CI) | EF5% (95% CI) | BEDROC ($\alpha=20$) | ROC-AUC | PR-AUC | Accuracy @ Top 1% |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **1. 1D Physicochemical (Control)** | `0.0x` `(0.0–0.0)` | `0.0x` `(0.0–0.8)` | `0.026` | `0.346` | `0.019` | `96.5%` |
| **2. 2D ECFP4 (vs. `ZM241385`)** | **`36.0x`** `(16.0–40.0)` | **`9.6x`** `(5.6–13.6)` | **`0.519`** | `0.725` | **`0.433`** | `98.3%` |
| **3. 3D Shape + USRCAT (`ETKDGv3` + `O3A`)** | `8.0x` `(0.0–20.0)` | `6.4x` `(3.2–10.4)` | `0.355` | `0.722` | `0.166` | `96.9%` |
| **4. 3D Docking (`Vina · 4EIY dry`)** | `0.0x` `(0.0–0.0)` | `0.8x` `(0.0–2.8)` | `0.063` | `0.584` | `0.033` | `96.5%` |
| **5. 3D Docking (`Vina · 4EIY hydrated`)** | `4.0x` `(0.0–12.0)` | `1.6x` `(0.0–4.0)` | `0.128` | `0.593` | `0.043` | `96.7%` |

- **Metric contrast:** `2D ECFP4` (`0.725`) and `3D Shape + USRCAT` (`0.722`) show nearly identical global `ROC-AUC`, but `2D ECFP4` recovers 4.5× more actives in the Top 1% (`EF1% = 36.0x` vs. `8.0x`). Predicting all compounds inactive yields `97.5%` accuracy.
- **Mechanistic interpretation of the A/B docking test:** while retaining the 10 ordered waters resolves cognate `ZM241385` redocking (`10/10` `< 2.0 Å`) and raises `EF1%` from `0.0x` to `4.0x` (`BEDROC` `0.063` $\rightarrow$ `0.128`), overall docking enrichment remains well below 2D and 3D ligand-based baselines. On the strict non-`ZM241385`-analogue subset (`n = 21` actives with `ECFP4 < 0.35`), `4EIY_hydrated` yields `EF1% = 0.0x` and lower global `ROC-AUC` (`0.571`) than `4EIY_dry` (`0.625`). Treating `ZM241385`'s crystallographic waters as non-displaceable steric atoms blocks `298 / 975` decoys ($\Delta G > 0\text{ kcal/mol}$) but simultaneously introduces steric clashes (`> 0 kcal/mol`) for `23 / 60` non-cognate $A_{2A}$ actives (`Tozadenant`, `Imaradenant`) that must displace first-shell waters upon binding, compounded by rigid single-crystal sidechains and Vina's approximate desolvation scoring.

![ROC and Early Enrichment Curves](reports/figures/benchmark_roc_and_enrichment.png)

![Accuracy and ROC-AUC vs Early Enrichment Metrics](reports/figures/metrics_accuracy_vs_enrichment.png)

---

### 4. Screening Score vs. Experimental Binding Affinity ($pK_i$)

Across **60 experimentally measured $K_i$ binders** spanning $pK_i \in [4.0,\; 11.0]$, neither `4EIY_dry` Vina ($r = +0.16$, $R^2 = 0.03$, $p = 0.22$), nor `4EIY_hydrated` Vina ($r = +0.12$, $R^2 = 0.01$, $p = 0.36$), nor 2D ECFP4 similarity ($r = +0.25$, $R^2 = 0.06$) shows meaningful linear correlation with experimental $pK_i$.

![Screening Score vs Experimental pKi](reports/figures/affinity_vs_score_scatter.png)

---

## Repository Structure

- [`src/adora2a_audit/`](src/adora2a_audit/)
  - [`structures.py`](src/adora2a_audit/structures.py) — GPCRdb + RCSB PDB audit, atomic pocket residue checks, and `inspect_pocket_waters` (`gemmi`).
  - [`curation.py`](src/adora2a_audit/curation.py) — `CHEMBL251` curation with provenance, RDKit standardization, and contextual PAINS/Brenk auditing.
  - [`decoys.py`](src/adora2a_audit/decoys.py) — 7D property-matched decoy generation and 1D physicochemical leakage gate.
  - [`metrics.py`](src/adora2a_audit/metrics.py) — tie-aware analytical implementations of `EF1%`, `EF5%`, `BEDROC` ($\alpha=20$), `ROC-AUC`, `PR-AUC`, and stratified bootstrap CIs.
  - [`benchmark.py`](src/adora2a_audit/benchmark.py) — dry/hydrated `4EIY` PDBQT preparation, 10-seed cognate redocking, 2D similarity + *Murcko scaffold split*, 3D shape (`USRCAT` + `CrippenO3A`), and parallel `--batch` AutoDock Vina scoring.
- [`scripts/run_audit_pipeline.py`](scripts/run_audit_pipeline.py) — deterministic CLI pipeline.
- [`notebooks/01_adora2a_retrospective_audit.ipynb`](notebooks/01_adora2a_retrospective_audit.ipynb) — executed walkthrough notebook.
- [`reports/structure_audit.md`](reports/structure_audit.md) — audit of all 188 `ADORA2A` structures and `4EIY` crystallographic waters.
- [`reports/executive_summary.md`](reports/executive_summary.md) — executive summary and pipeline recommendations for AI/data engineering teams.
- [`tests/`](tests/) — unit tests (`pytest`) covering curation, metrics, structure/water selection, hydrated receptor assembly, and symmetry-corrected RMSD.

---

## Reproduction

```bash
# 1. Install package in editable mode with dev dependencies
pip install -e ".[dev]"

# 2. Run the end-to-end audit pipeline (reuses cached files in data/cache/ and processed scores)
python scripts/run_audit_pipeline.py

# 3. Run unit tests and linter
pytest -v
ruff check .
```

## Data Sources & Licensing

- **GPCRdb:** receptor state and co-crystallized ligand annotations via the [GPCRdb REST API](https://gpcrdb.org/).
- **RCSB PDB:** crystallographic metadata and atomic coordinates subject to the [RCSB PDB Usage Policy](https://www.rcsb.org/pages/usage-policy).
- **EMBL-EBI ChEMBL:** bioactivity records and chemical structures licensed under [CC BY-SA 3.0](https://chembl.gitbook.io/chembl-interface-documentation/about).
