# Structural Audit of Human ADORA2A (`P29274` · `aa2ar_human`)

## 1. Purpose of the Audit

In structure-based virtual screening (SBVS) and 3D generative modeling on GPCRs, **naively selecting the lowest nominal resolution PDB entry (`ORDER BY resolution ASC LIMIT 1`) is a frequent and costly failure mode**:

1. **Conformational state mismatch (agonist vs. antagonist):** $G_s$-coupled activation of the human Adenosine $A_{2A}$ receptor (`ADORA2A`) triggers an ~14 Å outward swing of intracellular TM6, contraction of the deep ribose-binding subpocket ($\text{Thr88}^{3.36}$, $\text{Ser277}^{7.42}$, $\text{His278}^{7.43}$), and collapse of the conserved Class A allosteric $\text{Na}^+$ cavity at $\text{Asp52}^{2.50}$. Screening antagonists against active-state structures—or ribose agonists against inactive-state structures—introduces systematic geometric bias.
2. **Crystallographic engineering artefacts (`StaR2` vs. intracellular fusions):** many ultra-high-resolution $A_{2A}$ crystal structures rely on the thermostabilized **`A2A-StaR2`** construct (8–9 point mutations), which introduces **`T88A` ($3.36$)** and **`S277A` ($7.42$)** **directly inside the orthosteric pocket**, deleting two critical polar hydrogen-bond donors.
3. **Protonation state & conserved water network:** inside the orthosteric site, $\text{Glu169}^{\text{ECL2}}$ and $\text{Asn253}^{6.55}$ form directional hydrogen bonds (2.79–3.17 Å) with the heteroaromatic ligand core. Naive hydrogenation without pH 7.4 carboxylate deprotonation converts $\text{Glu169}$ into a neutral donor (`-COOH`), sterically clashing with `ZM241385`'s exocyclic amine and flipping the crystallographic pose by 180°.

---

## 2. Global Landscape Across All 188 Experimental Structures

By merging the live **GPCRdb** structure catalog (`/services/structure/protein/aa2ar_human/`) with polymer-entity and mutation metadata from **RCSB PDB** (`data/processed/adora2a_structures_audit.csv`), we audited **188 human experimental structures**:

| Conformational State (GPCRdb) | X-ray Diffraction | Cryo-EM | Electron Crystallography | Total | Resolved Allosteric $\text{Na}^+$ ($\text{Asp52}^{2.50}$) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Inactive** | 166 | 4 | 2 | **172** | **145 / 172 (84.3%)** *(100% at $\le 1.90\text{ Å}$)* |
| **Active** | 11 | 4 | 0 | **15** | **0 / 15 (0%)** — $\text{Na}^+$ pocket collapsed |
| **Intermediate** | 0 | 1 | 0 | **1** | **0 / 1 (0%)** (`9XQB`, inverse-agonist Fab `mAb16`) |
| **Total** | **177** | **9** | **2** | **188** | **145 / 188 (77.1%)** |

![ADORA2A Structural Landscape](figures/structure_audit_resolution_state.png)

> **Two quantitative takeaways across the full 188-entry catalog:**
> 1. **Pervasive pocket-mutation risk (`66.0%`):** **124 out of 188 structures (`66.0%`)** are thermostabilized `StaR2` constructs carrying direct risk of `T88A` and `S277A` pocket mutations, compared with only **39 structures (`20.7%`)** preserving a wild-type 7TM bundle stabilized solely via an intracellular `ICL3-BRIL` fusion.
> 2. **Allosteric $\text{Na}^+$ switch:** **0 out of 15 active structures (`0%`)** retain the $\text{Na}^+$ ion at $\text{Asp52}^{2.50}$, versus **145 out of 172 inactive structures (`84.3%`, and `100%` at $\le 1.90\text{ Å}$)**—directly reflecting allosteric $\text{Na}^+$ expulsion upon receptor activation.

---

## 3. Critical Comparison of High-Resolution Candidates

### 3.1 Candidates for **Antagonist / Inverse-Agonist** Screening (`Inactive` State)

| PDB | Resolution | Method | Orthosteric Ligand | Fusion / Complex | Receptor Chain Mutations | Atomic Pocket Verification (`structures.py`) | Decision |
| :--- | :---: | :--- | :--- | :--- | :---: | :--- | :--- |
| **`5NM4`** | **1.70 Å** | XFEL | `ZM241385` (`ZMA`) | `BRIL` (ICL3) | **9 (`StaR2`)** | Shifted sequential numbering + **`T88A` & `S277A` inside the pocket** | **Rejected** (mutated pocket) |
| **`5IU4`** | **1.72 Å** | X-ray | `ZM241385` (`ZMA`) | `BRIL` (ICL3) | **11 (`StaR2`)** | Verified in coordinates: `88 = ALA` (`T88A`) and `277 = ALA` (`S277A`) | **Rejected** (mutated pocket) |
| **`9H37`** | **1.72 Å** | X-ray | `Compound 9` (`A1IR0`) | `BRIL` (ICL3) | 0 in 7TM | Recent series-specific inverse agonist; narrower benchmark comparability | Secondary backup |
| **`4EIY`** | **1.80 Å** | X-ray | **`ZM241385` (`ZMA`)** | `BRIL` (ICL3, res. 1001–1106) | **0 in 7TM** *(3 inside BRIL only)* | **100% Wild-Type** (`D52`, `T88`, `F168`, `E169`, `W246`, `H250`, `N253`, `S277`, `H278`) + resolved $\text{Na}^+$ | **SELECTED (Inactive Reference)** |
| **`7IO5`** | **1.84 Å** | X-ray | `Theophylline` + fragment `FU12-45` | `BRIL` (`StaR2`) | **8 (`StaR2`)** | Co-soaked secondary fragments + `StaR2` pocket mutations | **Rejected** |
| **`3EML`** | **2.60 Å** | X-ray | `ZM241385` (`ZMA`) | `T4L` (ICL3) | 0 in 7TM | Lower resolution (2.60 Å); lacks resolved $\text{Na}^+$ ion and ordered deep waters | **Rejected** (superseded by `4EIY`) |

### 3.2 Candidates for **Agonist** Studies (`Active` State)

| PDB | Resolution | Method | Orthosteric Ligand | Signaling Partner / Fusion | Mutations | Structural Verification | Decision |
| :--- | :---: | :--- | :--- | :--- | :---: | :--- | :--- |
| **`9EE8`** | **2.63 Å** | Cryo-EM | **`Adenosine` (`ADN`)** | **Full heterotrimeric $G_s$ + `Nb35`** (no ICL3 fusion) | 2 (outside pocket) | **100% Wild-Type** orthosteric pocket; bona fide $G_s$-coupled physiological active state at 2.63 Å | **SELECTED (Active Reference)** |
| **`5G53`** | **3.40 Å** | X-ray | `NECA` (`NEC`) | `mini-Gs` | 1 (`N154Q` deglycosylation) | Wild-type pocket coupled to `mini-Gs`, but lower resolution (3.40 Å) than `9EE8` | Historical crystal benchmark |
| **`4UHR`** | **2.60 Å** | X-ray | `CGS 21680` (`NGI`) | None (`GL31` thermostabilized) | **5 (`GL31`)** | Lacks coupled G protein (*active-like intermediate* intracellularly) and carries 5 thermostabilizing mutations | **Rejected as primary active** |
| **`5WF5`** | **2.60 Å** | X-ray | `UK-432,097` (`UKA`) | `T4L` (ICL3), no G protein | 5 | Stabilized solely by bulky agonist + `T4L` fusion without $G_s$ coupling | **Rejected** |

---

## 4. Explicit Selection & Orthosteric Box Geometry for Phase 3

1. **Selected Inactive Reference (`4EIY`, 1.80 Å):**
   - Trading just **0.08–0.10 Å** of nominal resolution relative to `5NM4`/`5IU4` preserves wild-type `Thr88` and `Ser277` in the deep polar pocket.
   - The chimeric cytochrome $b_{562}\text{RIL}$ (`BRIL`, residues `1001–1106`) sits in intracellular loop `ICL3` **> 28 Å** away from the orthosteric centroid, allowing clean excision of `BRIL` without perturbing the 7TM bundle or extracellular disulfide-stabilized loops (`ECL2`/`ECL3`).
   - **Verified crystallographic recognition geometry of `ZM241385` in `4EIY`:**
     - `Glu169.OE2` $\cdots$ `ZMA.N15` (exocyclic amine): **2.79 Å**
     - `Asn253.OD1` $\cdots$ `ZMA.N15` (exocyclic amine): **2.79 Å**
     - `Asn253.ND2` $\cdots$ `ZMA.O25` (furan oxygen): **3.01 Å** and `ZMA.N17`: **3.17 Å**
     - Aromatic $\pi$-$\pi$ stacking against `Phe168` (`ECL2`).
   - **Orthosteric docking box (defined strictly from co-crystallized `ZMA` heavy atoms):** center `(-0.424, 8.527, 17.127)` Å, dimensions `20.0 × 20.0 × 20.0` Å.

2. **Selected Active Reference (`9EE8`, 2.63 Å):**
   - Captures the canonical $G_s$-coupled active conformation bound to endogenous **adenosine (`ADN`)** with a 100% wild-type pocket, establishing why ribofuranosyl agonists **must be excluded** from retrospective antagonist screens on inactive `4EIY`.
