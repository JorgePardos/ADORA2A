"""Unit tests for hydrated receptor preparation and symmetry-corrected in-place RMSD."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from adora2a_audit.benchmark import (
    _build_crystal_heavy_graph,
    compute_inplace_symmetry_rmsd,
    prepare_4eiy_hydrated_receptor,
    prepare_4eiy_receptor_and_reference_ligand,
)
from adora2a_audit.structures import inspect_pocket_waters

CACHE_DIR = Path(__file__).resolve().parents[1] / "data" / "cache"
PDB_4EIY = CACHE_DIR / "structures" / "pdb_00004eiy.pdb"


@pytest.mark.skipif(not PDB_4EIY.exists(), reason="Requires cached 4EIY PDB file")
def test_hydrated_receptor_preserves_exact_dry_protein_atoms_and_appends_10_oda_waters():
    waters = inspect_pocket_waters(PDB_4EIY)
    retained = [w for w in waters if w.retained]
    assert len(retained) == 10

    with tempfile.TemporaryDirectory() as tmpd:
        dry_pdbqt, _, _ = prepare_4eiy_receptor_and_reference_ligand(PDB_4EIY, tmpd)
        hyd_pdbqt = prepare_4eiy_hydrated_receptor(dry_pdbqt, retained, tmpd)

        dry_lines = dry_pdbqt.read_text(encoding="utf-8").splitlines()
        hyd_lines = hyd_pdbqt.read_text(encoding="utf-8").splitlines()

        # Exact protein lines must be 100% identical between dry and hydrated arms
        assert hyd_lines[: len(dry_lines)] == dry_lines
        appended = hyd_lines[len(dry_lines) :]
        # 10 retained waters * (1 OA + 2 HD) = 30 HETATM lines
        assert len(appended) == 30
        assert sum(line.endswith("OA") for line in appended) == 10
        assert sum(line.endswith("HD") for line in appended) == 20


@pytest.mark.skipif(not PDB_4EIY.exists(), reason="Requires cached 4EIY PDB file")
def test_inplace_symmetry_rmsd_is_zero_on_crystal_and_invariant_to_phenol_automorphism():
    cryst_graph, cryst_coords = _build_crystal_heavy_graph(PDB_4EIY, resname="ZMA")
    with tempfile.TemporaryDirectory() as tmpd:
        _, zma_3d_h, _ = prepare_4eiy_receptor_and_reference_ligand(PDB_4EIY, tmpd)
    zma_no_h = Chem.RemoveHs(zma_3d_h)

    # Exact crystal pose must give 0.0 A RMSD
    assert compute_inplace_symmetry_rmsd(zma_no_h, 0, cryst_graph, cryst_coords) == pytest.approx(
        0.0, abs=1e-4
    )

    # Rigid translation by +3.0 A along X must give exact 3.0 A in-place RMSD (not aligned away!)
    shifted = Chem.Mol(zma_no_h)
    conf = shifted.GetConformer(0)
    for i in range(shifted.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + 3.0, p.y, p.z))
    assert compute_inplace_symmetry_rmsd(shifted, 0, cryst_graph, cryst_coords) == pytest.approx(
        3.0, abs=1e-4
    )

    # Applying a non-identity graph automorphism (180-deg phenol ring flip) must still give 0.0 A!
    automorphisms = [
        m for m in zma_no_h.GetSubstructMatches(zma_no_h, uniquify=False)
        if list(m) != list(range(zma_no_h.GetNumAtoms()))
    ]
    assert len(automorphisms) >= 1
    perm = automorphisms[0]
    flipped = Chem.Mol(zma_no_h)
    fconf = flipped.GetConformer(0)
    orig_pos = np.array(zma_no_h.GetConformer(0).GetPositions())
    for dst_idx, src_idx in enumerate(perm):
        fconf.SetAtomPosition(dst_idx, orig_pos[src_idx].tolist())
    assert compute_inplace_symmetry_rmsd(flipped, 0, cryst_graph, cryst_coords) == pytest.approx(
        0.0, abs=1e-4
    )
