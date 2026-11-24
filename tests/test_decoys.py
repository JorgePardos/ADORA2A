"""Unit tests for physicochemical descriptor calculation and decoy bias / dissimilarity gates."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from adora2a_audit.decoys import DESCRIPTOR_COLS, compute_1d_descriptors_and_fp

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"


def test_compute_1d_descriptors_and_fp_on_caffeine():
    props, fp = compute_1d_descriptors_and_fp("Cn1cnc2c1c(=O)n(C)c(=O)n2C")
    assert props is not None and fp is not None
    assert set(props.keys()) == set(DESCRIPTOR_COLS)
    assert props["mw"] == pytest.approx(194.19, abs=0.2)
    assert props["formal_charge"] == 0
    assert fp.GetNumBits() == 2048


@pytest.mark.skipif(
    not (PROC_DIR / "benchmark_1000_scores.csv").exists(),
    reason="Requires generated benchmark_1000_scores.csv",
)
def test_generated_benchmark_enforces_dissimilarity_and_zero_leakage():
    df = pd.read_csv(PROC_DIR / "benchmark_1000_scores.csv")
    assert len(df) == 1000
    act = df[df["is_active"] == 1]
    dec = df[df["is_active"] == 0]
    assert len(act) == 25
    assert len(dec) == 975

    # Zero InChIKey connectivity overlap between actives and decoys
    act_conn = {ik[:14] for ik in act["inchikey"]}
    dec_conn = {ik[:14] for ik in dec["inchikey"]}
    assert act_conn.isdisjoint(dec_conn)

    # All decoys must have ECFP4 Tanimoto < 0.35 against co-crystallized active ZM241385
    assert float(dec["score_2d_zma"].max()) < 0.35
