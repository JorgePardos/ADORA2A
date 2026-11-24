"""Four-level retrospective virtual screening benchmark & potency correlation audit.

Compares on the exact same set of curated ADORA2A antagonists + property-matched decoys:
  1. 1D Physicochemical negative control (5-fold OOF RandomForest on 7 descriptors)
  2. 2D Topological similarity (ECFP4 Tanimoto vs 4EIY co-crystal ZMA + Murcko scaffold split)
  3. 3D Ligand shape/pharmacophore (ETKDGv3 conformers + USRCAT & O3A Shape Tanimoto vs ZMA)
  4. 3D Structure-based docking (AutoDock Vina v1.2.7 into wild-type 4EIY orthosteric pocket)

Also quantifies:
  - Cognate redocking accuracy & dehydration degeneracy on 4EIY (`ZM241385`)
  - Scaffold-hop degradation (Stratified Random CV vs Murcko Scaffold Group CV)
  - Docking score / similarity vs experimental `pKi` correlation (`score != affinity`)
"""

from __future__ import annotations

import contextlib
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import gemmi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import (
    AllChem,
    rdFingerprintGenerator,
    rdMolAlign,
    rdMolDescriptors,
    rdShapeHelpers,
)
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold

from adora2a_audit.metrics import (
    ScreeningMetrics,
    compute_roc_curve,
    enrichment_curve,
    evaluate_screening,
)

RDLogger.DisableLog("rdApp.*")

STD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


@dataclass(frozen=True)
class RedockingAudit:
    """Cognate redocking validation summary for 4EIY + ZM241385."""

    top_pose_kcal: float
    native_like_pose_rank: int
    native_like_pose_kcal: float
    delta_kcal: float
    furan_oxygen_distance_angstrom: float


def prepare_4eiy_receptor_and_reference_ligand(
    pdb_path: str | Path,
    out_dir: str | Path,
    zma_smiles: str = "Oc1ccc(cc1)CCNc1nc(N)n2c(n1)nc(n2)c1ccco1",
) -> tuple[Path, Chem.Mol, tuple[float, float, float]]:
    """Prepare pH 7.4 7TM receptor PDBQT from 4EIY (excluding ICL3-BRIL, waters, lipids)
    and extract the 3D crystallographic ZMA conformer with chemical bond orders.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    st = gemmi.read_structure(str(pdb_path))
    chain = st[0]["A"]

    clean_st = gemmi.Structure()
    clean_st.name = "4EIY_7TM"
    clean_model = gemmi.Model("1")
    clean_chain = gemmi.Chain("A")

    zma_pdb_lines: list[str] = []
    zma_coords: list[list[float]] = []
    het_idx = 1

    for res in chain:
        seqnum = res.seqid.num
        if res.name in STD_AA and -5 <= seqnum <= 320:
            new_res = gemmi.Residue()
            new_res.name = res.name
            new_res.seqid = res.seqid
            seen = set()
            for atom in res:
                if atom.name in seen:
                    continue
                seen.add(atom.name)
                new_res.add_atom(atom)
            clean_chain.add_residue(new_res)
        elif res.name == "ZMA":
            for atom in res:
                if atom.element.name != "H":
                    zma_coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
                    zma_pdb_lines.append(
                        f"HETATM{het_idx:5d} {atom.name:<4s} ZMA A2401    "
                        f"{atom.pos.x:8.3f}{atom.pos.y:8.3f}{atom.pos.z:8.3f}  1.00  0.00          {atom.element.name:>2s}"
                    )
                    het_idx += 1

    clean_model.add_chain(clean_chain)
    clean_st.add_model(clean_model)
    clean_pdb = out / "4eiy_7tm_clean.pdb"
    clean_st.write_pdb(str(clean_pdb))

    rec_mol = Chem.MolFromPDBFile(str(clean_pdb), removeHs=False, sanitize=True)
    rec_mol = Chem.AddHs(rec_mol, addCoords=True)
    AllChem.ComputeGasteigerCharges(rec_mol)
    conf = rec_mol.GetConformer()

    pdbqt_lines: list[str] = []
    for atom in rec_mol.GetAtoms():
        idx = atom.GetIdx()
        anum = atom.GetAtomicNum()
        info = atom.GetPDBResidueInfo()
        rname = info.GetResidueName().strip() if info else "UNK"
        if anum == 1:
            nbr = atom.GetNeighbors()[0]
            if nbr.GetAtomicNum() == 6:
                continue
            nbr_info = nbr.GetPDBResidueInfo()
            nbr_res = nbr_info.GetResidueName().strip() if nbr_info else ""
            nbr_aname = nbr_info.GetName().strip() if nbr_info else ""
            # Physiological pH 7.4: ASP and GLU carboxylates are anionic H-bond acceptors
            if (nbr_res == "ASP" and nbr_aname in {"OD1", "OD2"}) or (
                nbr_res == "GLU" and nbr_aname in {"OE1", "OE2"}
            ):
                continue
            ad_type = "HD"
        elif anum == 6:
            ad_type = "A" if atom.GetIsAromatic() else "C"
        elif anum == 7:
            has_polar_h = any(n.GetAtomicNum() == 1 for n in atom.GetNeighbors())
            ad_type = "N" if has_polar_h else "NA"
        elif anum == 8:
            ad_type = "OA"
        elif anum == 16:
            has_polar_h = any(n.GetAtomicNum() == 1 for n in atom.GetNeighbors())
            ad_type = "S" if has_polar_h else "SA"
        else:
            ad_type = atom.GetSymbol()

        aname = info.GetName() if info else f" {atom.GetSymbol():<3s}"
        rnum = info.GetResidueNumber() if info else 1
        chid = info.GetChainId() if info else "A"
        pos = conf.GetAtomPosition(idx)
        q = float(atom.GetProp("_GasteigerCharge")) if atom.HasProp("_GasteigerCharge") else 0.0
        if np.isnan(q) or np.isinf(q):
            q = 0.0
        pdbqt_lines.append(
            f"ATOM  {idx + 1:5d} {aname:4s} {rname:3s} {chid:1s}{rnum:4d}    "
            f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00    {q:+6.3f} {ad_type:<2s}"
        )

    rec_pdbqt = out / "4eiy_7tm_ph74.pdbqt"
    rec_pdbqt.write_text("\n".join(pdbqt_lines) + "\n", encoding="utf-8")

    # Build reference 3D ZMA RDKit Mol with bond orders and hydrogens for 3D shape/USRCAT
    raw_zma = Chem.MolFromPDBBlock("\n".join(zma_pdb_lines) + "\nEND\n", sanitize=False, removeHs=True)
    tmpl = Chem.MolFromSmiles(zma_smiles)
    zma_3d = AllChem.AssignBondOrdersFromTemplate(tmpl, raw_zma)
    Chem.SanitizeMol(zma_3d)
    zma_3d_h = Chem.AddHs(zma_3d, addCoords=True)

    c_arr = np.asarray(zma_coords, dtype=float).mean(axis=0)
    center = (round(float(c_arr[0]), 3), round(float(c_arr[1]), 3), round(float(c_arr[2]), 3))
    return rec_pdbqt, zma_3d_h, center


def audit_cognate_redocking(
    vina_exe: str | Path,
    receptor_pdbqt: str | Path,
    center: tuple[float, float, float],
    pdb_path: str | Path,
    out_dir: str | Path,
    zma_smiles: str = "Oc1ccc(cc1)CCNc1nc(N)n2c(n1)nc(n2)c1ccco1",
) -> RedockingAudit:
    """Redock ZM241385 from ab initio 3D SMILES embedding into dehydrated 4EIY."""
    out = Path(out_dir)
    mol = Chem.AddHs(Chem.MolFromSmiles(zma_smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol)

    prep = MoleculePreparation()
    pdbqt_str, _, _ = PDBQTWriterLegacy.write_string(prep.prepare(mol)[0])
    lig_in = out / "cognate_zma_in.pdbqt"
    lig_out = out / "cognate_zma_out.pdbqt"
    lig_in.write_text(pdbqt_str, encoding="utf-8")

    res = subprocess.run(
        [
            str(vina_exe),
            "--receptor", str(receptor_pdbqt),
            "--ligand", str(lig_in),
            "--center_x", f"{center[0]:.3f}",
            "--center_y", f"{center[1]:.3f}",
            "--center_z", f"{center[2]:.3f}",
            "--size_x", "20.0",
            "--size_y", "20.0",
            "--size_z", "20.0",
            "--exhaustiveness", "8",
            "--cpu", "4",
            "--seed", "42",
            "--out", str(lig_out),
        ],
        capture_output=True,
        text=True,
        check=True,
    )

    mode_energies: list[float] = []
    for line in res.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            mode_energies.append(float(parts[1]))

    # Extract crystallographic furan O25 and phenol O4 coordinates
    st = gemmi.read_structure(str(pdb_path))
    cryst_o25 = None
    cryst_o4 = None
    for r in st[0]["A"]:
        if r.name == "ZMA":
            for a in r:
                if a.name.strip() == "O25":
                    cryst_o25 = np.array([a.pos.x, a.pos.y, a.pos.z])
                elif a.name.strip() == "O4":
                    cryst_o4 = np.array([a.pos.x, a.pos.y, a.pos.z])

    pmol = PDBQTMolecule.from_file(str(lig_out), skip_typing=True)
    dock_mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pmol)[0])
    # Identify aromatic furan O vs phenolic OH in docked molecule
    furan_o_idx = next(
        a.GetIdx() for a in dock_mol.GetAtoms() if a.GetAtomicNum() == 8 and a.GetIsAromatic()
    )
    phenol_o_idx = next(
        a.GetIdx() for a in dock_mol.GetAtoms() if a.GetAtomicNum() == 8 and not a.GetIsAromatic()
    )

    best_rank = 1
    best_furan_d = 999.0
    for cid in range(dock_mol.GetNumConformers()):
        c = dock_mol.GetConformer(cid)
        pf = np.array(list(c.GetAtomPosition(furan_o_idx)))
        pp = np.array(list(c.GetAtomPosition(phenol_o_idx)))
        df = float(np.linalg.norm(pf - cryst_o25))
        dp = float(np.linalg.norm(pp - cryst_o4))
        if df < 2.5 and dp < 3.5:
            best_rank = cid + 1
            best_furan_d = df
            break
        if df < best_furan_d:
            best_rank = cid + 1
            best_furan_d = df

    top_e = mode_energies[0]
    nat_e = mode_energies[best_rank - 1]
    return RedockingAudit(
        top_pose_kcal=top_e,
        native_like_pose_rank=best_rank,
        native_like_pose_kcal=nat_e,
        delta_kcal=round(abs(nat_e - top_e), 3),
        furan_oxygen_distance_angstrom=round(best_furan_d, 2),
    )


def score_2d_ecfp4_similarity(
    smiles_list: list[str],
    query_smiles: str = "Oc1ccc(cc1)CCNc1nc(N)n2c(n1)nc(n2)c1ccco1",
) -> np.ndarray:
    """Compute Morgan radius-2 (ECFP4, 2048-bit) Tanimoto similarity against query ligand."""
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    q_fp = gen.GetFingerprint(Chem.MolFromSmiles(query_smiles))
    scores = []
    for smi in smiles_list:
        m = Chem.MolFromSmiles(smi)
        fp = gen.GetFingerprint(m)
        scores.append(float(DataStructs.TanimotoSimilarity(fp, q_fp)))
    return np.asarray(scores, dtype=float)


def evaluate_2d_scaffold_vs_random_split(
    smiles_list: list[str],
    murcko_list: list[str],
    y_true: np.ndarray,
    seed: int = 42,
    label_prefix: str = "2D ECFP4 ML",
) -> tuple[ScreeningMetrics, ScreeningMetrics]:
    """Compare supervised 2D ECFP4 LogisticRegression under Stratified Random CV vs
    Murcko Scaffold Group CV to quantify performance drop on unseen scaffolds.
    """
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    x_mat = np.vstack([gen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in smiles_list])
    yt = np.asarray(y_true, dtype=int)
    groups = np.asarray(murcko_list)

    # 1. Stratified Random 5-fold CV (optimistic when congeneric analogues exist)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    oof_rand = np.zeros(len(yt), dtype=float)
    for tr, te in skf.split(x_mat, yt):
        clf = LogisticRegression(C=1.0, max_iter=1000, random_state=seed)
        clf.fit(x_mat[tr], yt[tr])
        oof_rand[te] = clf.predict_proba(x_mat[te])[:, 1]

    # 2. Stratified Murcko Scaffold 5-fold Group CV (honest new-series generalization)
    sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
    oof_scaf = np.zeros(len(yt), dtype=float)
    for tr, te in sgkf.split(x_mat, yt, groups):
        clf = LogisticRegression(C=1.0, max_iter=1000, random_state=seed)
        clf.fit(x_mat[tr], yt[tr])
        oof_scaf[te] = clf.predict_proba(x_mat[te])[:, 1]

    return (
        evaluate_screening(yt, oof_rand, method=f"{label_prefix} (Random Split)"),
        evaluate_screening(yt, oof_scaf, method=f"{label_prefix} (Murcko Scaffold Split)"),
    )


def _shape_worker(args: tuple[int, str, bytes]) -> tuple[int, float, float]:
    """Compute USRCAT similarity and O3A 3D Shape Tanimoto against reference ZMA."""
    idx, smi, ref_pkl = args
    RDLogger.DisableLog("rdApp.*")
    ref_mol = Chem.Mol(ref_pkl)
    ref_usrcat = rdMolDescriptors.GetUSRCAT(ref_mol)
    ref_no_h = Chem.RemoveHs(ref_mol)
    ref_crippen = rdMolDescriptors._CalcCrippenContribs(ref_no_h)

    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return idx, 0.0, 0.0
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 1
    cids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=5, params=params))
    if not cids:
        params.useRandomCoords = True
        cids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=3, params=params))
    if not cids:
        return idx, 0.0, 0.0

    best_usr = 0.0
    best_shape = 0.0
    for cid in cids:
        with contextlib.suppress(RuntimeError, ValueError):
            AllChem.MMFFOptimizeMolecule(mol_h, confId=cid, maxIters=120)
        u_vec = rdMolDescriptors.GetUSRCAT(mol_h, confId=cid)
        u_sim = float(rdMolDescriptors.GetUSRScore(ref_usrcat, u_vec))
        best_usr = max(best_usr, u_sim)

    mol_no_h = Chem.RemoveHs(mol_h)
    prb_crippen = rdMolDescriptors._CalcCrippenContribs(mol_no_h)
    for cid in range(mol_no_h.GetNumConformers()):
        with contextlib.suppress(RuntimeError, ValueError):
            o3a = rdMolAlign.GetCrippenO3A(mol_no_h, ref_no_h, prb_crippen, ref_crippen, cid, 0)
            o3a.Align()
            shape_dist = rdShapeHelpers.ShapeTanimotoDist(mol_no_h, ref_no_h, confId1=cid, confId2=0)
            shape_sim = 1.0 - float(shape_dist)
            best_shape = max(best_shape, shape_sim)

    # Combined 3D ligand score (equal blend of 3D steric Shape Tanimoto + pharmacophoric USRCAT)
    return idx, best_usr, 0.5 * best_usr + 0.5 * best_shape


def score_3d_ligand_shape(
    smiles_list: list[str],
    ref_zma_mol_h: Chem.Mol,
    max_workers: int = 12,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate 3D conformers in parallel and compute USRCAT & 3D Shape Tanimoto vs ZMA."""
    ref_pkl = ref_zma_mol_h.ToBinary()
    tasks = [(i, smi, ref_pkl) for i, smi in enumerate(smiles_list)]
    usrcat_scores = np.zeros(len(smiles_list), dtype=float)
    combo_scores = np.zeros(len(smiles_list), dtype=float)

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        for idx, u_sc, c_sc in pool.map(_shape_worker, tasks, chunksize=16):
            usrcat_scores[idx] = u_sc
            combo_scores[idx] = c_sc
    return usrcat_scores, combo_scores


def _prep_pdbqt_worker(args: tuple[int, str]) -> tuple[int, str]:
    """Embed ligand from SMILES and prepare PDBQT string via Meeko."""
    idx, smi = args
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return idx, ""
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 1
    if AllChem.EmbedMolecule(mol_h, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol_h, params) != 0:
            return idx, ""
    with contextlib.suppress(RuntimeError, ValueError):
        AllChem.MMFFOptimizeMolecule(mol_h, maxIters=120)

    try:
        prep = MoleculePreparation()
        setups = prep.prepare(mol_h)
        if not setups:
            return idx, ""
        pdbqt_str, is_ok, _ = PDBQTWriterLegacy.write_string(setups[0])
        return (idx, pdbqt_str) if (is_ok and pdbqt_str.strip()) else (idx, "")
    except (RuntimeError, ValueError, KeyError, TypeError):
        return idx, ""


def _vina_batch_worker(
    args: tuple[int, list[tuple[int, str]], str, str, tuple[float, float, float], int, int],
) -> list[tuple[int, float]]:
    """Run AutoDock Vina `--batch` on a chunk of ligands so the receptor grid is computed once."""
    _chunk_id, items, vina_exe, rec_pdbqt, center, exhaustiveness, cpus = args
    results: list[tuple[int, float]] = []
    if not items:
        return results

    with tempfile.TemporaryDirectory() as tmpd:
        in_dir = Path(tmpd) / "in"
        out_dir = Path(tmpd) / "out"
        in_dir.mkdir()
        out_dir.mkdir()
        valid_indices: list[int] = []
        for idx, pdbqt_str in items:
            if not pdbqt_str:
                results.append((idx, 0.0))
                continue
            (in_dir / f"lig_{idx}.pdbqt").write_text(pdbqt_str, encoding="utf-8")
            valid_indices.append(idx)

        if valid_indices:
            subprocess.run(
                [
                    vina_exe,
                    "--receptor", rec_pdbqt,
                    "--batch", str(in_dir),
                    "--dir", str(out_dir),
                    "--center_x", f"{center[0]:.3f}",
                    "--center_y", f"{center[1]:.3f}",
                    "--center_z", f"{center[2]:.3f}",
                    "--size_x", "20.0",
                    "--size_y", "20.0",
                    "--size_z", "20.0",
                    "--exhaustiveness", str(exhaustiveness),
                    "--max_evals", "1500",
                    "--cpu", str(cpus),
                    "--num_modes", "1",
                    "--seed", "42",
                    "--verbosity", "0",
                ],
                capture_output=True,
                text=True,
                timeout=900,
                check=False,
            )
            for idx in valid_indices:
                out_file = out_dir / f"lig_{idx}_out.pdbqt"
                score_val = 0.0
                if out_file.exists():
                    for line in out_file.read_text(encoding="utf-8").splitlines():
                        if line.startswith("REMARK VINA RESULT:"):
                            parts = line.split()
                            if len(parts) >= 4:
                                score_val = float(parts[3])
                            break
                results.append((idx, score_val))
    return results


def score_autodock_vina(
    smiles_list: list[str],
    vina_exe: str | Path,
    receptor_pdbqt: str | Path,
    center: tuple[float, float, float],
    exhaustiveness: int = 4,
    workers: int = 8,
    cpus_per_worker: int = 2,
) -> np.ndarray:
    """Prepare 3D PDBQTs in parallel and dock via AutoDock Vina `--batch` chunks."""
    n = len(smiles_list)
    pdbqt_list = [""] * n
    with ProcessPoolExecutor(max_workers=12) as pool:
        for idx, pstr in pool.map(_prep_pdbqt_worker, enumerate(smiles_list), chunksize=25):
            pdbqt_list[idx] = pstr

    # Split into ~25-ligand chunks across `workers` parallel Vina batch processes
    n_chunks = min(max(n // 25, 1), n)
    indexed_items = list(enumerate(pdbqt_list))
    chunks = [indexed_items[i::n_chunks] for i in range(n_chunks)]
    tasks = [
        (cid, ch, str(vina_exe), str(receptor_pdbqt), center, exhaustiveness, cpus_per_worker)
        for cid, ch in enumerate(chunks)
    ]

    kcal = np.zeros(n, dtype=float)
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_vina_batch_worker, t): t[0] for t in tasks}
        for fut in as_completed(futures):
            batch_res = fut.result()
            for idx, val in batch_res:
                kcal[idx] = val
            done += len(batch_res)
            print(f"  [AutoDock Vina --batch] docked {done} / {n} molecules", flush=True)
    return kcal


def analyze_affinity_correlation(actives_scored_df: pd.DataFrame) -> dict[str, dict[str, float]]:
    """Quantify correlation between virtual screening scores and experimental pKi on actives."""
    pki = actives_scored_df["pki_median"].to_numpy(dtype=float)
    out: dict[str, dict[str, float]] = {}
    for col, label in [
        ("vina_neg_kcal", "AutoDock Vina (-ΔG, kcal/mol)"),
        ("score_2d_zma", "2D ECFP4 Tanimoto (vs ZM241385)"),
        ("score_3d_shape", "3D Shape + USRCAT (vs ZM241385)"),
    ]:
        vals = actives_scored_df[col].to_numpy(dtype=float)
        r, p_r = pearsonr(vals, pki)
        rho, p_rho = spearmanr(vals, pki)
        out[label] = {
            "pearson_r": round(float(r), 3),
            "pearson_r2": round(float(r**2), 3),
            "pearson_p": round(float(p_r), 4),
            "spearman_rho": round(float(rho), 3),
            "spearman_p": round(float(p_rho), 4),
        }
    return out


def plot_roc_and_early_enrichment(
    y_true: np.ndarray,
    method_scores: dict[str, np.ndarray],
    out_path: str | Path,
) -> None:
    """Plot full ROC curves alongside log-scaled early enrichment (cumulative recall vs top %)."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    colors = {
        "1D Fisicoquímico (Control)": "#7f7f7f",
        "2D ECFP4 (vs ZM241385)": "#1f77b4",
        "3D Forma + USRCAT (vs ZMA)": "#ff7f0e",
        "3D Docking (AutoDock Vina · 4EIY)": "#2ca02c",
    }

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.3), dpi=180)
    ax1.plot([0, 1], [0, 1], "--", color="#bbbbbb", lw=1.1, label="Azar (AUC = 0.50)")

    for name, scores in method_scores.items():
        m = evaluate_screening(y_true, scores, method=name)
        fpr, tpr = compute_roc_curve(y_true, scores)
        fx, fy = enrichment_curve(y_true, scores)
        col = colors.get(name, "#333333")

        ax1.plot(fpr, tpr, lw=2.1, color=col, label=f"{name} (AUC = {m.roc_auc:.2f})")
        ax2.plot(
            fx[1:] * 100.0,
            fy[1:] * 100.0,
            lw=2.1,
            color=col,
            label=f"{name} (EF1% = {m.ef1:.1f}x, BEDROC = {m.bedroc_20:.2f})",
        )

    # Random expectation on log-x plot
    x_grid = np.logspace(-1, 2, 200)
    ax2.plot(x_grid, x_grid, "--", color="#bbbbbb", lw=1.1, label="Azar (EF = 1.0x)")
    ax2.axvline(1.0, color="#d62728", linestyle=":", lw=1.0, alpha=0.75)
    ax2.axvline(5.0, color="#9467bd", linestyle=":", lw=1.0, alpha=0.75)
    ax2.text(1.08, 88, "Top 1%\n(EF1%)", fontsize=8, color="#d62728", fontweight="bold")
    ax2.text(5.4, 88, "Top 5%\n(EF5%)", fontsize=8, color="#9467bd", fontweight="bold")

    ax1.set_xlabel("Tasa de Falsos Positivos (FPR)", fontsize=10)
    ax1.set_ylabel("Tasa de Verdaderos Positivos (TPR)", fontsize=10)
    ax1.set_title("A. Curvas ROC globales (visión clásica ML)", fontsize=10.5, fontweight="bold")
    ax1.legend(loc="lower right", fontsize=8)
    ax1.grid(True, alpha=0.22)

    ax2.set_xscale("log")
    ax2.set_xlim(0.15, 100.0)
    ax2.set_ylim(0, 103)
    ax2.set_xlabel("% de quimioteca cribada (escala logarítmica)", fontsize=10)
    ax2.set_ylabel("% de activos recuperados", fontsize=10)
    ax2.set_title("B. Enriquecimiento temprano (lo que importa en laboratorio)", fontsize=10.5, fontweight="bold")
    ax2.legend(loc="lower right", fontsize=7.8)
    ax2.grid(True, which="both", alpha=0.22)

    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_metrics_contrast(
    metrics_list: list[ScreeningMetrics],
    ci_map: dict[str, dict[str, tuple[float, float]]],
    out_path: str | Path,
) -> None:
    """Contrast misleading metrics (Accuracy, ROC-AUC) vs honest early-enrichment metrics (BEDROC, EF1%)."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    names = [m.method for m in metrics_list]
    y_pos = np.arange(len(names))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.2), dpi=180)
    h = 0.25

    acc_vals = [m.accuracy_top1pct for m in metrics_list]
    auc_vals = [m.roc_auc for m in metrics_list]
    bed_vals = [m.bedroc_20 for m in metrics_list]
    pr_vals = [m.pr_auc for m in metrics_list]

    ax1.barh(y_pos + 1.5 * h, acc_vals, height=h, color="#cbd5e1", label="Accuracy @ Top 1% (engañosa)")
    ax1.barh(y_pos + 0.5 * h, auc_vals, height=h, color="#60a5fa", label="ROC-AUC (global)")
    ax1.barh(y_pos - 0.5 * h, bed_vals, height=h, color="#f59e0b", label="BEDROC (α=20, temprano)")
    ax1.barh(y_pos - 1.5 * h, pr_vals, height=h, color="#10b981", label="PR-AUC (sensible al desbalance)")

    naive_acc = metrics_list[0].naive_all_inactive_accuracy
    ax1.axvline(naive_acc, color="#64748b", linestyle="--", lw=1.0)
    ax1.text(naive_acc - 0.27, len(names) - 0.45, f"Predecir todo inactivo = {naive_acc*100:.1f}% acc!", fontsize=7.8, color="#475569")
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(names, fontsize=8.3)
    ax1.set_xlim(0, 1.08)
    ax1.set_xlabel("Valor de la métrica [0 – 1]", fontsize=10)
    ax1.set_title("A. Por qué Accuracy (~97.5%) y ROC-AUC ocultan fallos tempranos", fontsize=10.2, fontweight="bold")
    ax1.legend(loc="lower right", fontsize=7.8)
    ax1.grid(True, axis="x", alpha=0.22)

    ef1_vals = [m.ef1 for m in metrics_list]
    ef5_vals = [m.ef5 for m in metrics_list]
    ef1_err_lo = [m.ef1 - ci_map[m.method]["ef1"][0] for m in metrics_list]
    ef1_err_hi = [ci_map[m.method]["ef1"][1] - m.ef1 for m in metrics_list]

    ax2.barh(y_pos + 0.18, ef1_vals, xerr=[ef1_err_lo, ef1_err_hi], height=0.34, color="#ef4444", capsize=3, label="EF1% (con IC 95% bootstrap)")
    ax2.barh(y_pos - 0.18, ef5_vals, height=0.34, color="#8b5cf6", label="EF5%")
    ax2.axvline(1.0, color="#888888", linestyle="--", lw=1.0, label="Azar (EF = 1.0x)")

    for idx, (e1, e5, hi) in enumerate(zip(ef1_vals, ef5_vals, [ci_map[m.method]["ef1"][1] for m in metrics_list])):
        ax2.text(max(hi, e5) + 0.8, idx + 0.16, f"EF1%={e1:.1f}x", fontsize=8, va="center", fontweight="bold", color="#b91c1c")

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([""] * len(names))
    ax2.set_xlim(0, max(max(ci_map[m.method]["ef1"][1] for m in metrics_list) * 1.24, 14.0))
    ax2.set_xlabel("Factor de Enriquecimiento (veces sobre azar)", fontsize=10)
    ax2.set_title(f"B. Enriquecimiento temprano real (EF1% máx teórico = {metrics_list[0].ef1_max:.0f}x)", fontsize=10.2, fontweight="bold")
    ax2.legend(loc="lower right", fontsize=8)
    ax2.grid(True, axis="x", alpha=0.22)

    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_affinity_vs_scores(
    bench_df: pd.DataFrame,
    corr_stats: dict[str, dict[str, float]],
    out_path: str | Path,
) -> None:
    """Scatter of AutoDock Vina score and 2D similarity vs experimental pKi with annotated FP/FN."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    act = bench_df[bench_df["is_active"] == 1].copy()
    dec = bench_df[bench_df["is_active"] == 0].copy()

    top5_vina_cut = float(np.percentile(bench_df["vina_neg_kcal"], 95.0))
    top5_2d_cut = float(np.percentile(bench_df["score_2d_zma"], 95.0))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.4, 5.4), dpi=180)

    # Panel A: Vina (-kcal/mol) vs experimental pKi on actives + decoy distribution strip
    # Exclude rare rigid-box steric clashes (ΔG > 0 kcal/mol, i.e. -ΔG < 0) from visual y-range
    n_clashes = int((dec["vina_neg_kcal"] < 0.0).sum())
    dec_plot = dec[dec["vina_neg_kcal"] >= 0.0].copy()

    rng = np.random.default_rng(42)
    dec_x_jitter = rng.normal(loc=4.35, scale=0.11, size=len(dec_plot))
    ax1.scatter(
        dec_x_jitter,
        dec_plot["vina_neg_kcal"],
        s=15,
        alpha=0.22,
        color="#94a3b8",
        label=f"Decoys emparejados (n={len(dec_plot)}, +{n_clashes} choques estéricos)",
    )
    ax1.scatter(act["pki_median"], act["vina_neg_kcal"], s=58, alpha=0.88, color="#16a34a", edgecolors="white", linewidths=0.6, label=f"Activos curados Ki (n={len(act)})")
    ax1.axhline(top5_vina_cut, color="#dc2626", linestyle="--", lw=1.0, label=f"Corte Top 5% Vina ({-top5_vina_cut:.2f} kcal/mol)")

    # Stagger reference ligand labels cleanly
    ref_offsets = {
        "ZM241385": (6, 6),
        "Istradefylline": (6, -10),
        "Preladenant": (-58, 7),
        "Vipadenant": (6, 5),
        "Tozadenant": (-56, -11),
        "Imaradenant": (6, -10),
        "Caffeine": (6, 5),
        "Theophylline": (6, -10),
    }
    act["reference_name"] = act["reference_name"].fillna("").astype(str)
    for _, r in act[act["reference_name"] != ""].iterrows():
        rname = str(r["reference_name"])
        ax1.annotate(
            rname,
            xy=(r["pki_median"], r["vina_neg_kcal"]),
            xytext=ref_offsets.get(rname, (5, 4)),
            textcoords="offset points",
            fontsize=7.4,
            fontweight="bold",
            color="#0f5132",
        )

    fn_candidates = act[act["vina_neg_kcal"] < top5_vina_cut].sort_values("pki_median", ascending=False)
    if not fn_candidates.empty:
        fn = fn_candidates.iloc[0]
        lbl = fn["reference_name"] if fn["reference_name"] else str(fn["primary_chembl_id"])
        ax1.annotate(
            f"FN: {lbl}\n(pKi={fn['pki_median']:.2f}, ΔG={-fn['vina_neg_kcal']:.1f})",
            xy=(fn["pki_median"], fn["vina_neg_kcal"]),
            xytext=(-105, -26),
            textcoords="offset points",
            fontsize=7.3,
            color="#b91c1c",
            bbox={"boxstyle": "round,pad=0.2", "fc": "#fef2f2", "ec": "#ef4444"},
            arrowprops={"arrowstyle": "->", "color": "#ef4444", "lw": 0.9},
        )

    fp_decoy = dec_plot.sort_values("vina_neg_kcal", ascending=False).iloc[0]
    ax1.annotate(
        f"FP Decoy: {fp_decoy['primary_chembl_id']}\n(ΔG={-fp_decoy['vina_neg_kcal']:.1f} kcal/mol)",
        xy=(4.35, fp_decoy["vina_neg_kcal"]),
        xytext=(24, -4),
        textcoords="offset points",
        fontsize=7.3,
        color="#7c2d12",
        bbox={"boxstyle": "round,pad=0.2", "fc": "#fffbeb", "ec": "#f59e0b"},
        arrowprops={"arrowstyle": "->", "color": "#f59e0b", "lw": 0.9},
    )

    st_v = corr_stats["AutoDock Vina (-ΔG, kcal/mol)"]
    ax1.set_ylim(4.5, 13.0)
    ax1.set_title(
        f"A. Score AutoDock Vina vs. pKi experimental\n"
        f"(Activos: r = {st_v['pearson_r']:+.2f}, R² = {st_v['pearson_r2']:.2f}, ρ = {st_v['spearman_rho']:+.2f} → score ≠ afinidad)",
        fontsize=9.8,
        fontweight="bold",
    )
    ax1.set_xlabel("pKi experimental en ChEMBL (decoys en banda izquierda)", fontsize=9.5)
    ax1.set_ylabel("Afinidad predicha Vina −ΔG (kcal/mol, arriba = mejor)", fontsize=9.5)
    ax1.legend(loc="lower right", fontsize=7.6)
    ax1.grid(True, alpha=0.22)

    # Panel B: 2D ECFP4 similarity to ZM241385 vs experimental pKi
    dec_x_all = rng.normal(loc=4.35, scale=0.11, size=len(dec))
    ax2.scatter(dec_x_all, dec["score_2d_zma"], s=15, alpha=0.22, color="#94a3b8", label="Decoys (ECFP4 < 0.35 respecto a activos)")
    ax2.scatter(act["pki_median"], act["score_2d_zma"], s=58, alpha=0.88, color="#2563eb", edgecolors="white", linewidths=0.6, label="Activos curados Ki")
    ax2.axhline(top5_2d_cut, color="#dc2626", linestyle="--", lw=1.0, label=f"Corte Top 5% ECFP4 ({top5_2d_cut:.2f})")

    for _, r in act[act["reference_name"] != ""].iterrows():
        rname = str(r["reference_name"])
        ax2.annotate(
            rname,
            xy=(r["pki_median"], r["score_2d_zma"]),
            xytext=ref_offsets.get(rname, (5, 4)),
            textcoords="offset points",
            fontsize=7.4,
            fontweight="bold",
            color="#1e3a8a",
        )

    st_2d = corr_stats["2D ECFP4 Tanimoto (vs ZM241385)"]
    ax2.set_title(
        f"B. Similitud 2D ECFP4 (vs ZM241385) vs. pKi experimental\n"
        f"(Activos: r = {st_2d['pearson_r']:+.2f}, R² = {st_2d['pearson_r2']:.2f}, ρ = {st_2d['spearman_rho']:+.2f})",
        fontsize=9.8,
        fontweight="bold",
    )
    ax2.set_xlabel("pKi experimental en ChEMBL (decoys mostrados en banda izquierda)", fontsize=9.5)
    ax2.set_ylabel("Tanimoto ECFP4 frente a ZM241385", fontsize=9.5)
    ax2.legend(loc="upper left", fontsize=7.8)
    ax2.grid(True, alpha=0.22)

    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
