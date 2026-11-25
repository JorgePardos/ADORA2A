# Executive Summary · ADORA2A Retrospective Virtual Screening Audit

**Audience:** AI Engineering, Data Science, and Computational Chemistry Teams  
**Target:** Human Adenosine $A_{2A}$ Receptor (`ADORA2A`, UniProt `P29274`, ChEMBL `CHEMBL251`, Class A GPCR)  
**Scope & Objective:** Evaluate the upstream structural, chemical-curation, and benchmark-design controls required before deploying virtual screening or ML ranking models for $A_{2A}$ antagonists, including a controlled single-variable A/B test of crystallographic pocket hydration (`4EIY_dry` vs. `4EIY_hydrated`) in AutoDock Vina v1.2.7.

---

## 1. Explicit Upstream Quality Gates (*Go / No-Go*)

| Audit Gate | Quantitative *Go / No-Go* Criterion | Observed Result | Status |
| :--- | :--- | :--- | :---: |
| **1. 3D Pocket Integrity** | Zero engineered mutations in orthosteric/allosteric contact residues; conformational state matched to modality (`Inactive` for antagonists). | Excluded `5NM4` (1.70 Å) and `5IU4` (1.72 Å)—and **66.0% (`124/188`)** of $A_{2A}$ PDB entries—due to `StaR2` pocket mutations (`T88A`, `S277A`). Selected **`4EIY` (1.80 Å, wild-type 7TM pocket + allosteric $\text{Na}^+$)** for antagonists and **`9EE8` (2.63 Å, active + $G_s$)** for agonists. | **GO** |
| **2. Assay Purity & Provenance** | Retain `assay_chembl_id`, `confidence_score >= 8`, exact `'='` relation; separate equilibrium binding $K_i$ (`Ki_B`) from functional $IC_{50}$; exclude agonists from inactive screening. | From `12,105` raw records $\rightarrow$ `5,901` exact `Ki_B` measurements $\rightarrow$ **`4,972` unique parent structures (`InChIKey`)**, flagging `44` compounds with $> 1$ log-unit inter-assay spread and excluding `546` agonists. | **GO** |
| **3. Contextual Structural Alerts** | Audit PAINS/Brenk catalogs without automatic deletion of validated clinical series. | Blind PAINS/Brenk filtering would discard **3 co-crystallized Phase 2/3 clinical antagonists**: `Tozadenant` (PAINS `anil_di_alk_C`), `Vipadenant` (Brenk `aniline`), and `Imaradenant` (Brenk `2-halo_pyridine`). | **GO** |
| **4. 1D Decoy-Bias Control** | 5-fold CV 1D classifiers (`MW`, `cLogP`, `TPSA`, `HBD`, `HBA`, `RotBonds`, `Charge`) must yield $\text{ROC-AUC} \le 0.60$. | Random ChEMBL decoys show substantial 1D leakage (**AUC = `0.814` LogReg / `0.783` RF**, driven by higher `TPSA` in $A_{2A}$ ligands). Joint 7D matching + `ECFP4 < 0.35` reduces 1D AUC to **`0.484` (LogReg)** and **`0.346` (RF)**. Note that `0.346` in RF reflects slight inverted OOF signal and small-sample variance across 25 locally matched actives rather than exact `0.50` noise. | **GO** |

*Methodological note on receptor preparation:* Both `4EIY_dry` and `4EIY_hydrated` PDBQT files were generated at pH 7.4 (`Asp`/`Glu` deprotonated) via a custom `RDKit`/`gemmi` writer with Gasteiger partial charges; AutoDock Vina's `vina` scoring function does not use partial charges, relying only on AutoDock atom types and covalent `HD` connectivity (used to type the 10 retained waters—selected via $d_{\text{ZMA}} \le 4.5\text{ Å}$, $B < 30.0\text{ Å}^2$, occupancy $= 1.00$, $\ge 1$ H-bond $\le 3.5\text{ Å}$—as donor–acceptor `XS_TYPE_O_DA`).

---

## 2. Empirical Results: 10-Seed Redocking & 1,000-Compound Benchmark (2.5% Hits)

### 2.1 Cognate `ZM241385` Redocking Across 10 Independent Seeds
- **`4EIY_dry`:** Pose 1 RMSD = **`9.68 ± 0.02 Å`** (`0/10` seeds $< 2.0\text{ Å}$), despite sampling the native pose within the top 9 modes in **`9/10` seeds** (`Best RMSD = 1.34 ± 0.49 Å`). Across those 9 seeds, $\Delta E = E_{\text{native}} - E_{\text{non-native}} = \mathbf{+0.188 \pm 0.084\text{ kcal/mol}}$, demonstrating a reproducible scoring inversion rather than sampling failure.
- **`4EIY_hydrated` (`HOH 2517, 2518, 2521, 2524, 2525, 2527, 2531, 2541, 2583, 2585`):** Pose 1 RMSD = **`1.05 ± 0.15 Å`** (**`10/10` seeds $< 2.0\text{ Å}$**) and $\Delta E = \mathbf{-0.254 \pm 0.044\text{ kcal/mol}}$, resolving cognate pose inversion across all 10 seeds.

### 2.2 Retrospective Screening Benchmark (`25` Actives + `975` 7D-Matched Decoys · Naive Accuracy = `97.5%`)

| Level / Method | EF1% (95% CI) | EF5% (95% CI) | BEDROC ($\alpha=20$) | ROC-AUC | PR-AUC | Accuracy @ Top 1% | Engineering Takeaway |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **1. 1D Control (Physicochemical RF)** | `0.0x` `(0–0)` | `0.0x` `(0.0–0.8)` | `0.026` | `0.346` | `0.019` | `96.5%` | Passes `<0.60` bias gate (`0.346` shows slight local-matching inversion; `96.5%` accuracy is uninformative). |
| **2. 2D Baseline (ECFP4 vs. `ZM241385`)** | **`36.0x`** `(16–40)` | **`9.6x`** `(5.6–13.6)` | **`0.519`** | `0.725` | **`0.433`** | `98.3%` | Strongest early enrichment (`9/10` Top-1% hits), driven partly by congeneric analogues (`23.7x` without `ZMA` analogues). |
| **3. 3D Shape (`ETKDGv3` + `USRCAT` + `O3A`)** | `8.0x` `(0–20)` | `6.4x` `(3.2–10.4)` | `0.355` | `0.722` | `0.166` | `96.9%` | Matches 2D on global `ROC-AUC` (`0.722` vs. `0.725`), but yields 4.5× lower Top-1% recovery. |
| **4. 3D Docking (`Vina · 4EIY dry`)** | `0.0x` `(0–0)` | `0.8x` `(0.0–2.8)` | `0.063` | `0.584` | `0.033` | `96.5%` | Fails early enrichment once 1D MW/`TPSA` differences are removed (`AUC = 0.625` on non-`ZMA` analogues). |
| **5. 3D Docking (`Vina · 4EIY hydrated`)** | `4.0x` `(0–12)` | `1.6x` `(0.0–4.0)` | `0.128` | `0.593` | `0.043` | `96.7%` | Modest early gain (`BEDROC` `0.063` $\rightarrow$ `0.128`), but drops to `EF1% = 0.0x` (`AUC = 0.571`) on non-`ZMA` scaffolds. |

---

## 3. Data-Driven Mechanistic Diagnosis

1. **Why static pocket hydration solves cognate redocking yet yields only modest virtual-screening gains:**
   - Retaining `4EIY`'s 10 ordered pocket waters stabilizes `ZM241385`'s native pose (`10/10` seeds $< 2.0\text{ Å}$) and sterically excludes `298 / 975` (`30.6%`) matched decoys ($\Delta G > 0\text{ kcal/mol}$), doubling `BEDROC` from `0.063` to `0.128` and recovering `ZM241385`-compatible binders (`EF1% = 4.0x`).
   - However, because crystallographic waters in Vina are treated as **immutable steric obstacles** without water-displacement thermodynamics, those same 10 waters sterically clash ($\Delta G > 0\text{ kcal/mol}$) with **`23 / 60` (`38.3%`)** diverse $A_{2A}$ antagonists (including `Tozadenant` and `Imaradenant`) that displace subsets of first-shell waters upon binding. Consequently, on the 21 non-`ZM241385`-analogue actives, `4EIY_hydrated` drops back to `EF1% = 0.0x` (`ROC-AUC = 0.571` vs. `0.625` dry). Poor Vina enrichment against 7D-matched decoys therefore stems jointly from **single-conformation pocket rigidity**, **non-displaceable water representation**, and **Vina's approximate empirical contact/desolvation scoring**.
2. **Screening scores do not predict binding affinity ($pK_i$):**
   - Across **60 measured $K_i$ binders** ($pK_i \in [4.0,\; 11.0]$), linear variance explained is negligible for `4EIY_dry` Vina ($r = +0.161$, $R^2 = 0.026$, $p = 0.22$), `4EIY_hydrated` Vina ($r = +0.120$, $R^2 = 0.014$, $p = 0.36$), and 2D ECFP4 ($r = +0.246$, $R^2 = 0.060$).
3. **Scaffold splitting is mandatory for supervised 2D evaluation:**
   - On experimental `CHEMBL251` $K_i$ data ($n = 1,661$), switching from stratified random CV to **Murcko scaffold group CV** lowers `PR-AUC` from `0.917` to `0.880` and `BEDROC` from `0.947` to `0.914`.

---

## 4. Practical Recommendations for Production Pipelines

1. **Automate structural & pocket-mutation checks upstream:** reject `StaR2` pocket-mutated PDB entries (`T88A`, `S277A`) automatically, enforce pH 7.4 carboxylate protonation on `Glu169`/`Asp52`, and validate cognate pose recovery across $\ge 10$ seeds before screening.
2. **Use ligand-based 2D/3D triage ahead of displaceable-water / flexible-pocket rescoring:** prioritize fast 2D/3D ligand descriptors validated under Murcko scaffold splits for primary library reduction; when applying structure-based methods across diverse chemotypes, avoid freezing all ligand-specific first-shell waters rigidly and instead use displaceable-water protocols or physics-based rescoring.
3. **Enforce CI/CD benchmark gates:** report **EF1%**, **BEDROC ($\alpha=20$)**, and **PR-AUC** with bootstrap CIs, never *accuracy*, and run an automated 1D physicochemical adversarial check (`LogisticRegression` + `RandomForest` $\text{ROC-AUC} \le 0.60$) on every benchmark dataset.
