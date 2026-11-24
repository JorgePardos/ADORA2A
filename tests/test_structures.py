"""Unit tests for ADORA2A structural audit and pocket residue verification."""

from __future__ import annotations

from pathlib import Path

import pytest

from adora2a_audit.structures import build_structure_audit_table, inspect_pocket_coordinates

CACHE_DIR = Path(__file__).resolve().parents[1] / "data" / "cache"


@pytest.mark.skipif(
    not (CACHE_DIR / "gpcrdb_aa2ar_human.json").exists(),
    reason="Requires cached GPCRdb JSON in data/cache/",
)
def test_structure_catalog_separates_states_and_na_allostery():
    df = build_structure_audit_table(CACHE_DIR)
    assert len(df) >= 180
    assert {"Inactive", "Active", "Intermediate"}.issubset(set(df["state"]))

    # Active structures systematically lack the allosteric Na+ ion (pocket collapsed)
    active_df = df[df["state"] == "Active"]
    assert active_df["has_allosteric_na"].sum() == 0

    # 4EIY is present as Inactive 1.80 A WT-7TM + ICL3-BRIL with Na+
    r_4eiy = df[df["pdb_code"] == "4EIY"].iloc[0]
    assert r_4eiy["state"] == "Inactive"
    assert r_4eiy["resolution"] == pytest.approx(1.8)
    assert bool(r_4eiy["has_allosteric_na"]) is True
    assert r_4eiy["construct_class"] == "WT-7TM + ICL3-BRIL"


@pytest.mark.skipif(
    not (CACHE_DIR / "structures" / "pdb_00004eiy.cif").exists(),
    reason="Requires cached 4EIY/5IU4 coordinate files",
)
def test_pocket_coordinate_inspection_catches_star2_mutations_in_5iu4_vs_wt_4eiy():
    p_4eiy = inspect_pocket_coordinates(
        CACHE_DIR / "structures" / "pdb_00004eiy.cif", pdb_code="4EIY", ligand_resname="ZMA"
    )
    assert p_4eiy.is_wt_orthosteric_pocket is True
    assert p_4eiy.mutated_pocket_residues == {}
    assert p_4eiy.has_resolved_na_ion is True
    assert p_4eiy.ligand_centroid_xyz is not None

    p_5iu4 = inspect_pocket_coordinates(
        CACHE_DIR / "structures" / "pdb_00005iu4.cif", pdb_code="5IU4", ligand_resname="ZMA"
    )
    assert p_5iu4.is_wt_orthosteric_pocket is False
    # StaR2 mutates T88A and S277A directly inside the orthosteric pocket!
    assert p_5iu4.mutated_pocket_residues[88] == ("THR", "ALA")
    assert p_5iu4.mutated_pocket_residues[277] == ("SER", "ALA")
