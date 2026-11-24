"""Unit tests for virtual screening metrics (EF1%, EF5%, BEDROC, ROC-AUC, PR-AUC, accuracy trap)."""

from __future__ import annotations

import numpy as np
import pytest

from adora2a_audit.metrics import (
    bedroc_score,
    enrichment_factor,
    evaluate_screening,
    max_enrichment_factor,
)


def test_perfect_ranking_reaches_theoretical_maxima():
    # 20 actives in 1000 compounds (2% hit rate)
    y_true = np.array([1] * 20 + [0] * 980)
    y_score = np.linspace(10.0, 0.0, 1000)

    m = evaluate_screening(y_true, y_score, method="perfect")
    assert m.roc_auc == pytest.approx(1.0)
    assert m.pr_auc == pytest.approx(1.0)
    assert m.bedroc_20 == pytest.approx(1.0, abs=1e-4)
    # Top 1% = 10 compounds -> all 10 are actives -> observed rate 1.0 / 0.02 = 50.0
    assert m.ef1 == pytest.approx(50.0)
    assert m.ef1 == pytest.approx(m.ef1_max)
    # Top 5% = 50 compounds -> contains all 20 actives -> (20/50) / (20/1000) = 20.0
    assert m.ef5 == pytest.approx(20.0)
    assert m.ef5 == pytest.approx(m.ef5_max)


def test_worst_ranking_gives_zero_early_enrichment():
    y_true = np.array([1] * 20 + [0] * 980)
    y_score = np.linspace(0.0, 10.0, 1000)  # actives ranked dead last

    m = evaluate_screening(y_true, y_score, method="worst")
    assert m.roc_auc == pytest.approx(0.0)
    assert m.ef1 == pytest.approx(0.0)
    assert m.ef5 == pytest.approx(0.0)
    assert m.bedroc_20 < 1e-4


def test_constant_tie_scores_equal_random_expectation():
    # All molecules receive identical score -> tie-aware EF and BEDROC must be order-invariant
    y_true_front = np.array([1] * 25 + [0] * 975)
    y_true_back = np.array([0] * 975 + [1] * 25)
    y_score = np.ones(1000) * 0.5

    assert enrichment_factor(y_true_front, y_score, fraction=0.01) == pytest.approx(1.0)
    assert enrichment_factor(y_true_back, y_score, fraction=0.01) == pytest.approx(1.0)
    assert enrichment_factor(y_true_front, y_score, fraction=0.05) == pytest.approx(1.0)
    assert bedroc_score(y_true_front, y_score, alpha=20.0) == pytest.approx(
        bedroc_score(y_true_back, y_score, alpha=20.0), abs=1e-12
    )
    assert evaluate_screening(y_true_front, y_score).roc_auc == pytest.approx(0.5)


def test_late_recovery_trap_decent_roc_auc_zero_ef1_and_low_bedroc():
    """Illustrates why ROC-AUC alone misleads in virtual screening:
    If all 20 actives sit at ranks 201..220 out of 1000 (top 20%-22%),
    ROC-AUC is ~0.79 ('good' in generic ML), yet EF1% = 0, EF5% = 0, and BEDROC20 is tiny!
    """
    y_true = np.zeros(1000, dtype=int)
    y_true[200:220] = 1
    y_score = np.linspace(1.0, 0.0, 1000)

    m = evaluate_screening(y_true, y_score, method="late_recovery")
    assert m.roc_auc > 0.78
    assert m.ef1 == 0.0
    assert m.ef5 == 0.0
    assert m.bedroc_20 < 0.05
    # Meanwhile naive all-inactive accuracy is 98.0%!
    assert m.naive_all_inactive_accuracy == pytest.approx(0.98)


def test_validation_rejects_bad_inputs():
    with pytest.raises(ValueError, match="Length mismatch"):
        enrichment_factor([1, 0], [0.5])
    with pytest.raises(ValueError, match="at least one active"):
        enrichment_factor([0, 0, 0], [0.1, 0.2, 0.3])
    with pytest.raises(ValueError, match="NaN"):
        bedroc_score([1, 0], [float("nan"), 0.2])
    assert max_enrichment_factor(1000, 20, 0.01) == pytest.approx(50.0)
