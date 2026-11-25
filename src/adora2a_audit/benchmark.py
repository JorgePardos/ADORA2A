"""Retrospective virtual screening benchmark, A/B pocket hydration audit, & pKi correlation.

Compares on the exact same set of curated ADORA2A antagonists + 7D-matched decoys:
  1. 1D Physicochemical negative control (5-fold OOF RandomForest & LogisticRegression)
  2. 2D Topological similarity (ECFP4 Tanimoto vs 4EIY co-crystal ZMA + Murcko scaffold split)
  3. 3D Ligand shape/pharmacophore (ETKDGv3 conformers + USRCAT & O3A Shape Tanimoto vs ZMA)
  4. 3D Structure-based docking — Dry arm (AutoDock Vina v1.2.7 on dehydrated WT 4EIY)
  5. 3D Structure-based docking — Hydrated arm (AutoDock Vina v1.2.7 on 4EIY retaining the
     10 ordered first-shell pocket/bridging waters selected by `inspect_pocket_waters`)

Also quantifies:
  - 10-seed cognate redocking reproducibility (Pose 1 & best-pose symmetry-corrected RMSD,
    < 2.0 A recovery rate, and ΔE between native and non-native poses on Dry vs Hydrated 4EIY)
  - Scaffold-hop degradation (Stratified Random CV vs Murcko Scaffold Group CV)
  - Docking score / similarity vs experimental `pKi` correlation across 60 measured binders
"""

from __future__ import annotations

import contextlib
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
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
from adora2a_audit.structures import PocketWater

RDLogger.DisableLog("rdApp.*")

STD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

DEFAULT_REDOCK_SEEDS = (101, 202, 303, 404, 505, 606, 707, 808, 909, 1010)


@dataclass(frozen=True)
class MultiSeedRedockingAudit:
    """10-seed cognate redocking validation summary for ZM241385 on a receptor arm."""

    receptor_arm: str
    n_seeds: int
    pose1_rmsd_mean: float
    pose1_rmsd_std: float
    best_rmsd_mean: float
    best_rmsd_std: float
    pose1_lt_2a_count: int
    any_pose_lt_2a_count: int
    pose1_kcal_mean: float
    pose1_kcal_std: float
    delta_e_native_minus_nonnative_mean: float | None
    delta_e_native_minus_nonnative_std: float | None
    n_seeds_with_both_states: int
    per_seed_pose1_rmsd: tuple[float, ...]
    per_seed_best_rmsd: tuple[float, ...]
    per_seed_delta_e: tuple[float, ...]

    def to_dict(self) -> dict:
        return asdict(self)


def prepare_4eiy_receptor_and_reference_ligand(
    pdb_path: str | Path,
    out_dir: str | Path,
    zma_smiles: str = "Oc1ccc(cc1)CCNc1nc(N)n2c(n1)nc(n2)c1ccco1",
) -> tuple[Path, Chem.Mol, tuple[float, float, float]]:
    """Prepare pH 7.4 dry 7TM receptor PDBQT from 4EIY via custom RDKit/gemmi writer
    with Gasteiger charges (note: Vina's scoring function ignores partial charges and uses
    AutoDock atom types + HD donor connectivity) and extract 3D crystallographic ZMA.
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


def prepare_4eiy_hydrated_receptor(
    dry_receptor_pdbqt: str | Path,
    retained_waters: list[PocketWater],
    out_dir: str | Path,
) -> Path:
    """Create `4eiy_7tm_ph74_hydrated.pdbqt` by appending retained crystallographic waters.

    Uses the exact same protein PDBQT lines as `dry_receptor_pdbqt` and appends each
    retained water oxygen (`OA`) with two covalent `HD` atoms (`0.957 A` O-H bond length)
    so AutoDock Vina's internal C++ receptor parser promotes the water oxygen from
    acceptor-only (`XS_TYPE_O_A`) to donor-acceptor (`XS_TYPE_O_DA`).
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    lines = Path(dry_receptor_pdbqt).read_text(encoding="utf-8").splitlines()
    serial = len(lines) + 1

    for w in retained_waters:
        if not w.retained:
            continue
        wx, wy, wz = w.xyz
        lines.append(
            f"HETATM{serial:5d}  O   HOH {w.chain_id:1s}{w.seqid:4d}    "
            f"{wx:8.3f}{wy:8.3f}{wz:8.3f}  1.00{w.b_factor:6.2f}    -0.834 OA"
        )
        serial += 1
        lines.append(
            f"HETATM{serial:5d}  H1  HOH {w.chain_id:1s}{w.seqid:4d}    "
            f"{wx + 0.957:8.3f}{wy:8.3f}{wz:8.3f}  1.00  0.00    +0.417 HD"
        )
        serial += 1
        lines.append(
            f"HETATM{serial:5d}  H2  HOH {w.chain_id:1s}{w.seqid:4d}    "
            f"{wx - 0.240:8.3f}{wy + 0.927:8.3f}{wz:8.3f}  1.00  0.00    +0.417 HD"
        )
        serial += 1

    hyd_path = out / "4eiy_7tm_ph74_hydrated.pdbqt"
    hyd_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return hyd_path


def _build_crystal_heavy_graph(pdb_path: str | Path, resname: str = "ZMA") -> tuple[Chem.Mol, np.ndarray]:
    """Build element-typed single-bond topological graph + heavy-atom coordinate array from PDB."""
    st = gemmi.read_structure(str(pdb_path))
    lig_res = next(r for r in st[0]["A"] if r.name == resname)
    rw = Chem.RWMol()
    coords: list[list[float]] = []
    for atom in lig_res:
        if atom.element.name != "H":
            rd_atom = Chem.Atom(atom.element.name)
            rd_atom.SetNoImplicit(True)
            rw.AddAtom(rd_atom)
            coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
    c_arr = np.asarray(coords, dtype=float)
    for i in range(len(c_arr)):
        for j in range(i + 1, len(c_arr)):
            if float(np.linalg.norm(c_arr[i] - c_arr[j])) < 1.75:
                rw.AddBond(i, j, Chem.BondType.SINGLE)
    return rw.GetMol(), c_arr


def compute_inplace_symmetry_rmsd(
    dock_mol_no_h: Chem.Mol,
    conf_id: int,
    cryst_graph: Chem.Mol,
    cryst_coords: np.ndarray,
) -> float:
    """Compute exact in-place symmetry-corrected heavy-atom RMSD (without rigid-body alignment)."""
    dm = Chem.RWMol(dock_mol_no_h)
    for bond in dm.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in dm.GetAtoms():
        atom.SetIsAromatic(False)

    matches = cryst_graph.GetSubstructMatches(dm.GetMol(), uniquify=False)
    if not matches:
        raise RuntimeError("Graph isomorphism failed between docked pose and crystal reference.")
    d_pos = dock_mol_no_h.GetConformer(conf_id).GetPositions()
    return min(
        float(np.sqrt(((cryst_coords[list(m)] - d_pos) ** 2).sum(axis=1).mean()))
        for m in matches
    )


def _redock_seed_worker(
    args: tuple[int, str, str, str, tuple[float, float, float], bytes, np.ndarray],
) -> tuple[int, float, float, float, float | None]:
    """Run one cognate ZMA redocking seed and return (seed, pose1_rmsd, best_rmsd, pose1_kcal, delta_e)."""
    seed, vina_exe, rec_pdbqt, lig_pdbqt_str, center, cryst_pkl, cryst_coords = args
    RDLogger.DisableLog("rdApp.*")
    cryst_graph = Chem.Mol(cryst_pkl)

    with tempfile.TemporaryDirectory() as tmpd:
        lig_in = Path(tmpd) / "zma_in.pdbqt"
        lig_out = Path(tmpd) / "zma_out.pdbqt"
        lig_in.write_text(lig_pdbqt_str, encoding="utf-8")
        res = subprocess.run(
            [
                vina_exe,
                "--receptor", rec_pdbqt,
                "--ligand", str(lig_in),
                "--center_x", f"{center[0]:.3f}",
                "--center_y", f"{center[1]:.3f}",
                "--center_z", f"{center[2]:.3f}",
                "--size_x", "20.0",
                "--size_y", "20.0",
                "--size_z", "20.0",
                "--exhaustiveness", "8",
                "--cpu", "2",
                "--num_modes", "9",
                "--seed", str(seed),
                "--out", str(lig_out),
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        energies = [
            float(line.split()[1])
            for line in res.stdout.splitlines()
            if len(line.split()) >= 2 and line.split()[0].isdigit()
        ]
        pmol = PDBQTMolecule.from_file(str(lig_out), skip_typing=True)
        dmol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pmol)[0])
        rmsds = [
            compute_inplace_symmetry_rmsd(dmol, cid, cryst_graph, cryst_coords)
            for cid in range(dmol.GetNumConformers())
        ]

    nat_es = [e for r, e in zip(rmsds, energies, strict=False) if r < 2.0]
    non_es = [e for r, e in zip(rmsds, energies, strict=False) if r >= 2.0]
    delta_e = (min(nat_es) - min(non_es)) if (nat_es and non_es) else None
    return seed, rmsds[0], min(rmsds), energies[0], delta_e


def audit_cognate_redocking_multiseed(
    vina_exe: str | Path,
    receptor_pdbqt: str | Path,
    center: tuple[float, float, float],
    pdb_path: str | Path,
    receptor_arm: str,
    seeds: tuple[int, ...] = DEFAULT_REDOCK_SEEDS,
    zma_smiles: str = "Oc1ccc(cc1)CCNc1nc(N)n2c(n1)nc(n2)c1ccco1",
    max_workers: int = 5,
) -> MultiSeedRedockingAudit:
    """Evaluate cognate ZM241385 redocking across 10 independent random seeds."""
    mol = Chem.AddHs(Chem.MolFromSmiles(zma_smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol)

    prep = MoleculePreparation()
    lig_pdbqt_str, _, _ = PDBQTWriterLegacy.write_string(prep.prepare(mol)[0])
    cryst_graph, cryst_coords = _build_crystal_heavy_graph(pdb_path, resname="ZMA")
    cryst_pkl = cryst_graph.ToBinary()

    tasks = [
        (s, str(vina_exe), str(receptor_pdbqt), lig_pdbqt_str, center, cryst_pkl, cryst_coords)
        for s in seeds
    ]
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        rows = list(pool.map(_redock_seed_worker, tasks))

    rows.sort(key=lambda x: x[0])
    p1_rmsds = [r[1] for r in rows]
    best_rmsds = [r[2] for r in rows]
    p1_kcals = [r[3] for r in rows]
    del_es = [r[4] for r in rows if r[4] is not None]

    return MultiSeedRedockingAudit(
        receptor_arm=receptor_arm,
        n_seeds=len(seeds),
        pose1_rmsd_mean=round(float(np.mean(p1_rmsds)), 3),
        pose1_rmsd_std=round(float(np.std(p1_rmsds)), 3),
        best_rmsd_mean=round(float(np.mean(best_rmsds)), 3),
        best_rmsd_std=round(float(np.std(best_rmsds)), 3),
        pose1_lt_2a_count=sum(r < 2.0 for r in p1_rmsds),
        any_pose_lt_2a_count=sum(r < 2.0 for r in best_rmsds),
        pose1_kcal_mean=round(float(np.mean(p1_kcals)), 3),
        pose1_kcal_std=round(float(np.std(p1_kcals)), 3),
        delta_e_native_minus_nonnative_mean=round(float(np.mean(del_es)), 3) if del_es else None,
        delta_e_native_minus_nonnative_std=round(float(np.std(del_es)), 3) if del_es else None,
        n_seeds_with_both_states=len(del_es),
        per_seed_pose1_rmsd=tuple(round(x, 2) for x in p1_rmsds),
        per_seed_best_rmsd=tuple(round(x, 2) for x in best_rmsds),
        per_seed_delta_e=tuple(round(x, 3) for x in del_es),
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

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    oof_rand = np.zeros(len(yt), dtype=float)
    for tr, te in skf.split(x_mat, yt):
        clf = LogisticRegression(C=1.0, max_iter=1000, random_state=seed)
        clf.fit(x_mat[tr], yt[tr])
        oof_rand[te] = clf.predict_proba(x_mat[te])[:, 1]

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
    _, items, vina_exe, rec_pdbqt, center, exhaustiveness, cpus = args
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
    tag: str = "4EIY",
) -> np.ndarray:
    """Prepare 3D PDBQTs in parallel and dock via AutoDock Vina `--batch` chunks."""
    n = len(smiles_list)
    pdbqt_list = [""] * n
    with ProcessPoolExecutor(max_workers=12) as pool:
        for idx, pstr in pool.map(_prep_pdbqt_worker, enumerate(smiles_list), chunksize=25):
            pdbqt_list[idx] = pstr

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
            print(f"  [AutoDock Vina --batch · {tag}] docked {done} / {n} molecules", flush=True)
    return kcal


def analyze_affinity_correlation(actives_scored_df: pd.DataFrame) -> dict[str, dict[str, float]]:
    """Quantify correlation between virtual screening scores and experimental pKi on actives."""
    pki = actives_scored_df["pki_median"].to_numpy(dtype=float)
    candidates = [
        ("vina_neg_kcal", "AutoDock Vina · 4EIY dry (-ΔG, kcal/mol)"),
        ("vina_hyd_neg_kcal", "AutoDock Vina · 4EIY hydrated (-ΔG, kcal/mol)"),
        ("score_2d_zma", "2D ECFP4 Tanimoto (vs ZM241385)"),
        ("score_3d_shape", "3D Shape + USRCAT (vs ZM241385)"),
    ]
    out: dict[str, dict[str, float]] = {}
    for col, label in candidates:
        if col not in actives_scored_df.columns:
            continue
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
        "1D Physicochemical (Control)": "#7f7f7f",
        "2D ECFP4 (vs ZM241385)": "#1f77b4",
        "3D Shape + USRCAT (vs ZMA)": "#ff7f0e",
        "3D Docking (Vina · 4EIY dry)": "#2ca02c",
        "3D Docking (Vina · 4EIY hydrated)": "#0d9488",
    }
    linestyles = {
        "3D Docking (Vina · 4EIY dry)": "--",
        "3D Docking (Vina · 4EIY hydrated)": "-",
    }

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.4), dpi=180)
    ax1.plot([0, 1], [0, 1], ":", color="#bbbbbb", lw=1.1, label="Random (AUC = 0.50)")

    for name, scores in method_scores.items():
        m = evaluate_screening(y_true, scores, method=name)
        fpr, tpr = compute_roc_curve(y_true, scores)
        fx, fy = enrichment_curve(y_true, scores)
        col = colors.get(name, "#333333")
        ls = linestyles.get(name, "-")

        ax1.plot(fpr, tpr, lw=2.1, linestyle=ls, color=col, label=f"{name} (AUC = {m.roc_auc:.2f})")
        ax2.plot(
            fx[1:] * 100.0,
            fy[1:] * 100.0,
            lw=2.1,
            linestyle=ls,
            color=col,
            label=f"{name} (EF1% = {m.ef1:.1f}x, BEDROC = {m.bedroc_20:.2f})",
        )

    x_grid = np.logspace(-1, 2, 200)
    ax2.plot(x_grid, x_grid, ":", color="#bbbbbb", lw=1.1, label="Random (EF = 1.0x)")
    ax2.axvline(1.0, color="#d62728", linestyle=":", lw=1.0, alpha=0.75)
    ax2.axvline(5.0, color="#9467bd", linestyle=":", lw=1.0, alpha=0.75)
    ax2.text(1.06, 54, "Top 1%\n(EF1%)", fontsize=8, color="#d62728", fontweight="bold")
    ax2.text(5.3, 54, "Top 5%\n(EF5%)", fontsize=8, color="#9467bd", fontweight="bold")

    ax1.set_xlabel("False Positive Rate (FPR)", fontsize=10)
    ax1.set_ylabel("True Positive Rate (TPR)", fontsize=10)
    ax1.set_title("A. Global ROC Curves (Classical ML View)", fontsize=10.5, fontweight="bold")
    ax1.legend(loc="lower right", fontsize=7.6)
    ax1.grid(True, alpha=0.22)

    ax2.set_xscale("log")
    ax2.set_xlim(0.15, 100.0)
    ax2.set_ylim(0, 103)
    ax2.set_xlabel("% of Library Screened (log scale)", fontsize=10)
    ax2.set_ylabel("% of Actives Recovered", fontsize=10)
    ax2.set_title("B. Early Enrichment (Top 1% / Top 5% Recovery)", fontsize=10.5, fontweight="bold")
    ax2.legend(loc="upper left", fontsize=7.4, framealpha=0.92)
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

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.8, 5.6), dpi=180)
    h = 0.20

    acc_vals = [m.accuracy_top1pct for m in metrics_list]
    auc_vals = [m.roc_auc for m in metrics_list]
    bed_vals = [m.bedroc_20 for m in metrics_list]
    pr_vals = [m.pr_auc for m in metrics_list]

    ax1.barh(y_pos + 1.5 * h, acc_vals, height=h, color="#cbd5e1", label="Accuracy @ Top 1% (misleading)")
    ax1.barh(y_pos + 0.5 * h, auc_vals, height=h, color="#60a5fa", label="ROC-AUC (global)")
    ax1.barh(y_pos - 0.5 * h, bed_vals, height=h, color="#f59e0b", label="BEDROC (α=20, early)")
    ax1.barh(y_pos - 1.5 * h, pr_vals, height=h, color="#10b981", label="PR-AUC (imbalance-aware)")

    naive_acc = metrics_list[0].naive_all_inactive_accuracy
    ax1.axvline(naive_acc, color="#64748b", linestyle="--", lw=1.0)
    ax1.text(
        naive_acc - 0.34,
        len(names) - 0.62,
        f"Predict all inactive = {naive_acc * 100:.1f}% acc",
        fontsize=7.5,
        color="#475569",
    )
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(names, fontsize=8.2)
    ax1.set_xlim(0, 1.08)
    ax1.set_xlabel("Metric Value [0 – 1]", fontsize=10)
    ax1.set_title("A. Why Accuracy (~97.5%) and ROC-AUC Hide Early Failures", fontsize=10.2, fontweight="bold")
    ax1.legend(loc="lower right", fontsize=7.6)
    ax1.grid(True, axis="x", alpha=0.22)

    ef1_vals = [m.ef1 for m in metrics_list]
    ef5_vals = [m.ef5 for m in metrics_list]
    ef1_err_lo = [m.ef1 - ci_map[m.method]["ef1"][0] for m in metrics_list]
    ef1_err_hi = [ci_map[m.method]["ef1"][1] - m.ef1 for m in metrics_list]

    ax2.barh(
        y_pos + 0.17,
        ef1_vals,
        xerr=[ef1_err_lo, ef1_err_hi],
        height=0.32,
        color="#ef4444",
        capsize=3,
        label="EF1% (95% bootstrap CI)",
    )
    ax2.barh(y_pos - 0.17, ef5_vals, height=0.32, color="#8b5cf6", label="EF5%")
    ax2.axvline(1.0, color="#888888", linestyle="--", lw=1.0, label="Random (EF = 1.0x)")

    ci_uppers = [ci_map[m.method]["ef1"][1] for m in metrics_list]
    for idx, (e1, e5, hi) in enumerate(zip(ef1_vals, ef5_vals, ci_uppers, strict=False)):
        ax2.text(
            max(hi, e5) + 0.8,
            idx + 0.16,
            f"EF1%={e1:.1f}x",
            fontsize=8,
            va="center",
            fontweight="bold",
            color="#b91c1c",
        )

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([""] * len(names))
    ax2.set_xlim(0, max(max(ci_uppers) * 1.24, 14.0))
    ax2.set_xlabel("Enrichment Factor (fold over random)", fontsize=10)
    ax2.set_title(
        f"B. Early Enrichment Comparison (Theoretical EF1% max = {metrics_list[0].ef1_max:.0f}x)",
        fontsize=10.2,
        fontweight="bold",
    )
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
    """3-panel scatter comparing Vina 4EIY dry, Vina 4EIY hydrated, and 2D ECFP4 vs experimental pKi."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    act = bench_df[bench_df["is_active"] == 1].copy()
    dec = bench_df[bench_df["is_active"] == 0].copy()
    act["reference_name"] = act["reference_name"].fillna("").astype(str)

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2), dpi=180)
    rng = np.random.default_rng(42)

    ref_offsets = {
        "ZM241385": (5, 5),
        "Istradefylline": (5, -9),
        "Preladenant": (-52, 6),
        "Vipadenant": (5, 4),
        "Tozadenant": (-50, -10),
        "Imaradenant": (5, -9),
        "Caffeine": (5, 4),
        "Theophylline": (5, -9),
    }

    vina_panels = [
        (axes[0], "vina_neg_kcal", "AutoDock Vina · 4EIY dry (-ΔG, kcal/mol)", "#16a34a", "A. Vina (4EIY dry) vs. pKi"),
        (axes[1], "vina_hyd_neg_kcal", "AutoDock Vina · 4EIY hydrated (-ΔG, kcal/mol)", "#0d9488", "B. Vina (4EIY hydrated) vs. pKi"),
    ]

    for ax, col, stat_key, dot_col, title_prefix in vina_panels:
        top5_cut = float(np.percentile(bench_df[col], 95.0))
        n_clashes = int((dec[col] < 0.0).sum())
        dec_plot = dec[dec[col] >= 0.0].copy()
        act_plot = act[act[col] >= 0.0].copy()
        act_clashes = int((act[col] < 0.0).sum())

        dec_x = rng.normal(loc=4.35, scale=0.11, size=len(dec_plot))
        clash_note = f"+{n_clashes} clashes" + (f", {act_clashes} active clash" if act_clashes else "")
        ax.scatter(
            dec_x,
            dec_plot[col],
            s=13,
            alpha=0.22,
            color="#94a3b8",
            label=f"Matched decoys (n={len(dec_plot)}, {clash_note})",
        )
        ax.scatter(
            act_plot["pki_median"],
            act_plot[col],
            s=52,
            alpha=0.88,
            color=dot_col,
            edgecolors="white",
            linewidths=0.6,
            label=f"Curated Ki actives (n={len(act_plot)})",
        )
        ax.axhline(top5_cut, color="#dc2626", linestyle="--", lw=1.0, label=f"Top 5% cutoff ({-top5_cut:.2f} kcal/mol)")

        for _, r in act_plot[act_plot["reference_name"] != ""].iterrows():
            rname = str(r["reference_name"])
            ax.annotate(
                rname,
                xy=(r["pki_median"], r[col]),
                xytext=ref_offsets.get(rname, (5, 4)),
                textcoords="offset points",
                fontsize=7.0,
                fontweight="bold",
                color="#0f5132",
            )

        fn_cands = act_plot[act_plot[col] < top5_cut].sort_values("pki_median", ascending=False)
        if not fn_cands.empty:
            fn = fn_cands.iloc[0]
            lbl = fn["reference_name"] if fn["reference_name"] else str(fn["primary_chembl_id"])
            fn_offset = (-112, 16) if float(fn[col]) < 5.0 else (-98, -24)
            ax.annotate(
                f"FN: {lbl}\n(pKi={fn['pki_median']:.2f}, ΔG={-fn[col]:.1f})",
                xy=(fn["pki_median"], fn[col]),
                xytext=fn_offset,
                textcoords="offset points",
                fontsize=6.9,
                color="#b91c1c",
                bbox={"boxstyle": "round,pad=0.2", "fc": "#fef2f2", "ec": "#ef4444"},
                arrowprops={"arrowstyle": "->", "color": "#ef4444", "lw": 0.85},
            )

        if not dec_plot.empty:
            fp_dec = dec_plot.sort_values(col, ascending=False).iloc[0]
            ax.annotate(
                f"FP: {fp_dec['primary_chembl_id']}\n(ΔG={-fp_dec[col]:.1f} kcal/mol)",
                xy=(4.35, fp_dec[col]),
                xytext=(20, -6),
                textcoords="offset points",
                fontsize=6.9,
                color="#7c2d12",
                bbox={"boxstyle": "round,pad=0.2", "fc": "#fffbeb", "ec": "#f59e0b"},
                arrowprops={"arrowstyle": "->", "color": "#f59e0b", "lw": 0.85},
            )

        st_v = corr_stats[stat_key]
        y_hi = max(13.0, float(bench_df[col].max()) + 0.8)
        ax.set_ylim(3.5, min(y_hi, 14.5))
        ax.set_title(
            f"{title_prefix}\n(r = {st_v['pearson_r']:+.2f}, R² = {st_v['pearson_r2']:.2f}, ρ = {st_v['spearman_rho']:+.2f})",
            fontsize=9.4,
            fontweight="bold",
        )
        ax.set_xlabel("Experimental ChEMBL pKi", fontsize=9.0)
        ax.set_ylabel("Predicted −ΔG (kcal/mol, higher = stronger)", fontsize=9.0)
        ax.legend(loc="upper right", fontsize=7.0, framealpha=0.92)
        ax.grid(True, alpha=0.22)

    # Panel C: 2D ECFP4 similarity to ZM241385 vs experimental pKi
    ax3 = axes[2]
    top5_2d_cut = float(np.percentile(bench_df["score_2d_zma"], 95.0))
    dec_x_all = rng.normal(loc=4.35, scale=0.11, size=len(dec))
    ax3.scatter(dec_x_all, dec["score_2d_zma"], s=13, alpha=0.22, color="#94a3b8", label="Decoys (ECFP4 < 0.35)")
    ax3.scatter(act["pki_median"], act["score_2d_zma"], s=52, alpha=0.88, color="#2563eb", edgecolors="white", linewidths=0.6, label="Curated Ki actives")
    ax3.axhline(top5_2d_cut, color="#dc2626", linestyle="--", lw=1.0, label=f"Top 5% cutoff ({top5_2d_cut:.2f})")

    for _, r in act[act["reference_name"] != ""].iterrows():
        rname = str(r["reference_name"])
        ax3.annotate(
            rname,
            xy=(r["pki_median"], r["score_2d_zma"]),
            xytext=ref_offsets.get(rname, (5, 4)),
            textcoords="offset points",
            fontsize=7.0,
            fontweight="bold",
            color="#1e3a8a",
        )

    st_2d = corr_stats["2D ECFP4 Tanimoto (vs ZM241385)"]
    ax3.set_title(
        f"C. 2D ECFP4 (vs ZM241385) vs. pKi\n(r = {st_2d['pearson_r']:+.2f}, R² = {st_2d['pearson_r2']:.2f}, ρ = {st_2d['spearman_rho']:+.2f})",
        fontsize=9.4,
        fontweight="bold",
    )
    ax3.set_xlabel("Experimental ChEMBL pKi", fontsize=9.0)
    ax3.set_ylabel("ECFP4 Tanimoto to ZM241385", fontsize=9.0)
    ax3.legend(loc="upper left", fontsize=7.2)
    ax3.grid(True, alpha=0.22)

    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
