# Executive Summary · ADORA2A Retrospective Virtual Screening Audit

**Target Audience:** AI Engineering, Data Science, and Computational Chemistry Teams  
**Biological Target:** Human Adenosine $A_{2A}$ Receptor (`ADORA2A`, UniProt `P29274`, ChEMBL `CHEMBL251`, Class A GPCR)  
**Core Engineering Question:** *Before scaling a virtual screening pipeline or training an ML potency model across multi-million-compound libraries for $A_{2A}$ antagonists (immuno-oncology / Parkinson's disease), which structural and data-curation gates prevent benchmark leakage, and which methods deliver genuine early enrichment once trivial shortcuts are removed?*

---

## 1. Explicit Quality Gates (*Go / No-Go*)

| Audit Gate | Quantitative *Go / No-Go* Criterion | Empirical Finding | Status |
| :--- | :--- | :--- | :---: |
| **1. 3D Pocket Integrity** | Zero engineered mutations in orthosteric/allosteric contact residues; conformational state matched to ligand modality (`Inactive` for antagonists). | Rejected `5NM4` (1.70 Å) and `5IU4` (1.72 Å)—along with **66.0% (`124/188`)** of all $A_{2A}$ PDB entries—due to `StaR2` mutations (`T88A`, `S277A`) **inside the pocket**. Selected **`4EIY` (1.80 Å, 100% wild-type pocket + allosteric $\text{Na}^+$)** for antagonists and **`9EE8` (2.63 Å, active + $G_s$)** for agonists. | **GO** |
| **2. Assay Purity & Provenance** | Retain `assay_chembl_id`, `confidence_score >= 8`, exact `'='` relation; separate equilibrium binding $K_i$ (`Ki_B`) from functional $IC_{50}$; exclude ribose agonists from inactive screening. | From `12,105` raw records $\rightarrow$ `5,901` exact `Ki_B` measurements $\rightarrow$ **`4,972` unique parent structures (`InChIKey`)**, flagging `44` compounds with $> 1$ log-unit inter-assay discrepancy and excluding `546` agonists. | **GO** |
| **3. Contextual Structural Alerts** | Audit PAINS/Brenk alerts without blind deletion of clinically validated chemotypes. | Blind PAINS/Brenk filtering would have silently deleted **3 co-crystallized Phase 2/3 clinical antagonists**: `Tozadenant` (PAINS `anil_di_alk_C`), `Vipadenant` (Brenk `aniline`), and `Imaradenant` (Brenk `2-halo_pyridine`). | **GO** |
| **4. 1D Decoy-Bias Control** | A 1D physicochemical classifier (`MW`, `cLogP`, `TPSA`, `HBD`, `HBA`, `RotBonds`, `Charge`) must yield cross-validated $\text{ROC-AUC} \le 0.60$. | Unmatched ChEMBL decoys suffer severe 1D leakage (**AUC = `0.81` LogReg / `0.78` RF**, driven by elevated `TPSA` in $A_{2A}$ heterocycles). Joint 7D matching + `ECFP4 < 0.35` reduces 1D AUC to **`0.48` / `0.35`**. | **GO** |

---

## 2. Honest Retrospective Benchmark Results (1,000 Compounds · 2.5% Hit Rate)

Evaluated on the exact same 1D-debiased benchmark (`25` scaffold-diverse potent antagonists + `975` 7D-matched decoys; baseline hit rate = `2.5%`, naive "predict all inactive" accuracy = **`97.5%`**, theoretical $\text{EF1\%}_{\max} = 40.0\times$):

| Level / Method | EF1% (95% CI) | EF5% | BEDROC ($\alpha=20$) | ROC-AUC | PR-AUC | Accuracy @ Top 1% | Key Takeaway for Engineering |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **1. 1D Control (Physicochemical RF)** | `0.0x` `(0–0)` | `0.0x` | `0.026` | `0.346` | `0.019` | `96.5%` | Confirms absence of 1D property leakage; note the deceptive `96.5%` accuracy. |
| **2. 2D Baseline (ECFP4 vs. `ZM241385`)** | **`36.0x`** `(16–40)` | **`9.6x`** | **`0.519`** | `0.725` | **`0.433`** | `98.3%` | Recovers 9/10 hits in Top 1%, but relies on topological proximity to query (`23.7x` once close analogues are removed). |
| **3. 3D Shape (`ETKDGv3` + `USRCAT` + `O3A`)** | `8.0x` `(0–20)` | `6.4x` | `0.355` | `0.722` | `0.166` | `96.9%` | **Identical global ROC-AUC (`0.722` vs. `0.725`) to 2D**, yet **4.5× lower Top-1% enrichment**: textbook proof that ROC-AUC is insufficient. |
| **4. 3D Docking (`AutoDock Vina` · `4EIY`)** | **`0.0x`** `(0–0)` | `0.8x` | `0.063` | `0.584` | `0.033` | `96.5%` | **Clear negative result:** once MW/lipophilicity shortcuts are removed, rigid water-free docking fails in the Top 1%. |

---

## 3. What Failed and Why (Mechanistic Diagnostics)

1. **Dehydrated rigid-receptor docking (`AutoDock Vina`) underperforms 2D similarity and 3D ligand shape:**
   - In our cognate `ZM241385` redocking audit on water-stripped `4EIY`, a head-to-tail flipped pose (`-9.08 kcal/mol`) and the native crystal orientation (rank 3, `-8.93 kcal/mol`, furan oxygen within `0.41 Å` of X-ray) differ by only **`0.15 kcal/mol`**—well inside Vina's scoring noise. Stripping conserved bridging waters (`HOH 2503–2505`) removes polar desolvation penalties in the deep cavity, allowing matched polyaromatic decoys (`-10` to `-12 kcal/mol`) to crowd the Top 1% (false positives) while rigid `ZM241385`-adapted sidechains penalize bulkier clinical antagonists (false negatives, e.g., `CHEMBL2419139`, $pK_i = 11.0$, $\Delta G = -6.6\text{ kcal/mol}$).
2. **Virtual screening score $\neq$ experimental binding affinity ($pK_i$):**
   - Across **60 experimentally measured $K_i$ ligands** spanning > 5 log units ($pK_i \in [4.0,\; 11.0]$), AutoDock Vina explains only **2.6% of variance** ($r = +0.161$, $R^2 = 0.026$, $p = 0.22$) and 2D ECFP4 explains **6.0%** ($r = +0.246$, $R^2 = 0.060$). Fast screening scores rank enriched hit lists; they are never drop-in surrogates for $\Delta G_{\text{bind}}$ (which requires FEP/ABFE or biophysical assays).
3. **Scaffold-hop degradation under *Murcko scaffold splitting*:**
   - On experimental `CHEMBL251` $K_i$ data ($n = 1,661$), switching from stratified random CV to **Murcko scaffold group CV** drops `PR-AUC` from `0.917` to `0.880` and `BEDROC` from `0.947` to `0.914`; similarly, removing close `ZM241385` analogues (`ECFP4 >= 0.35`) cuts unsupervised 2D `EF1%` from `36.0x` to `23.7x`.

---

## 4. Actionable Recommendations for Production Pipelines

1. **Automated pocket-mutation & protonation gate:** programmatically cross-check GPCRdb state + atomic residue identities before receptor preparation (blocking `StaR2` `T88A`/`S277A` constructs and enforcing pH 7.4 anionic carboxylates on `Glu169`/`Asp52`).
2. **Hierarchical funnel architecture (Fast 2D/3D $\rightarrow$ Hydrated Ensemble Docking $\rightarrow$ Physics Rescoring):** triage libraries with Murcko-validated 2D models and 3D shape/pharmacophore descriptors (`USRCAT`/`O3A`), then dock candidates into a **multi-crystal inactive ensemble** (`4EIY`, `5OLZ`, `6GT3`) **retaining ordered deep-pocket waters** and enforcing H-bond pharmacophore constraints on `Asn253^{6.55}` and `Phe168^{ECL2}`.
3. **Metric governance in ML CI/CD:** ban *accuracy*, require **EF1%**, **BEDROC ($\alpha=20$)**, and **PR-AUC** with bootstrap confidence intervals, and enforce an automated **1D adversarial decoy gate ($\text{ROC-AUC} < 0.60$)** on every benchmark update.
