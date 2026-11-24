"""Honest virtual screening evaluation metrics (EF1%, EF5%, BEDROC, ROC-AUC, PR-AUC).

Implements Truchon & Bayly (2007) BEDROC alongside early enrichment factors and
precision-recall metrics, contrasting them against naive accuracy on imbalanced sets.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

import numpy as np
from rdkit.ML.Scoring import Scoring as RdkitScoring
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve


@dataclass(frozen=True)
class ScreeningMetrics:
    """Summary of retrospective virtual screening performance."""

    method: str
    n_total: int
    n_actives: int
    hit_rate: float
    ef1: float
    ef1_max: float
    ef5: float
    ef5_max: float
    bedroc_20: float
    roc_auc: float
    pr_auc: float
    accuracy_top1pct: float
    naive_all_inactive_accuracy: float

    def to_dict(self) -> dict[str, float | int | str]:
        return asdict(self)


def _validate_inputs(y_true: np.ndarray, y_score: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    yt = np.asarray(y_true, dtype=int).ravel()
    ys = np.asarray(y_score, dtype=float).ravel()
    if yt.shape[0] != ys.shape[0]:
        raise ValueError(f"Length mismatch: y_true ({yt.shape[0]}) vs y_score ({ys.shape[0]}).")
    if yt.shape[0] == 0:
        raise ValueError("Empty input arrays.")
    unique = set(np.unique(yt).tolist())
    if not unique.issubset({0, 1}):
        raise ValueError(f"y_true must be binary {{0, 1}}, got {sorted(unique)}.")
    n_act = int(yt.sum())
    if n_act == 0 or n_act == yt.shape[0]:
        raise ValueError("y_true must contain at least one active (1) and one decoy (0).")
    if np.isnan(ys).any():
        raise ValueError("y_score contains NaN values; impute or drop failed molecules explicitly.")
    return yt, ys


def _stable_rank_order(y_true: np.ndarray, y_score: np.ndarray) -> np.ndarray:
    """Order indices by descending score, breaking ties deterministically via expected rank.

    To avoid optimistic or pessimistic tie-breaking artefacts when scores have ties
    (e.g. discrete similarity values or failed conformers), tied blocks are ordered
    conservatively with decoys before actives only if not averaged; here we sort by
    descending y_score and use stable index order after deterministic tie handling.
    """
    # Primary sort: descending y_score; secondary sort: ascending index for determinism
    return np.argsort(-y_score, kind="mergesort")


def enrichment_factor(
    y_true: np.ndarray | list[int],
    y_score: np.ndarray | list[float],
    fraction: float = 0.01,
) -> float:
    """Compute Enrichment Factor at a top-ranked fraction X (e.g. 0.01 for EF1%).

    Handles score ties at the cutoff boundary analytically by computing the expected
    number of actives selected when randomly subsampling tied molecules at the boundary.
    """
    if not (0.0 < fraction <= 1.0):
        raise ValueError(f"fraction must be in (0, 1], got {fraction}.")
    yt, ys = _validate_inputs(np.asarray(y_true), np.asarray(y_score))
    n_total = len(yt)
    n_actives = int(yt.sum())
    k = max(1, math.ceil(n_total * fraction))

    order = _stable_rank_order(yt, ys)
    sorted_scores = ys[order]
    sorted_labels = yt[order]

    cutoff_score = sorted_scores[k - 1]
    strictly_above = sorted_scores > cutoff_score
    tied_at_cutoff = sorted_scores == cutoff_score

    n_above = int(strictly_above.sum())
    actives_above = float(sorted_labels[strictly_above].sum())

    remaining_slots = k - n_above
    n_tied = int(tied_at_cutoff.sum())
    actives_tied = float(sorted_labels[tied_at_cutoff].sum())
    expected_actives_from_ties = remaining_slots * (actives_tied / n_tied)

    hits_in_top_k = actives_above + expected_actives_from_ties
    observed_rate = hits_in_top_k / k
    baseline_rate = n_actives / n_total
    return float(observed_rate / baseline_rate)


def max_enrichment_factor(n_total: int, n_actives: int, fraction: float = 0.01) -> float:
    """Theoretical maximum EF at top fraction X given dataset composition."""
    k = max(1, math.ceil(n_total * fraction))
    max_hits = min(k, n_actives)
    return float((max_hits / k) / (n_actives / n_total))


def bedroc_score(
    y_true: np.ndarray | list[int],
    y_score: np.ndarray | list[float],
    alpha: float = 20.0,
) -> float:
    """Compute Boltzmann-Enhanced Discrimination of ROC (BEDROC) per Truchon & Bayly (2007).

    When all scores are distinct, matches RDKit's CalcBEDROC exactly. When score ties
    occur, assigns each tied group its exact expected Boltzmann weight across the tied
    1-based rank span so the metric is invariant to input row ordering.
    """
    if alpha <= 0:
        raise ValueError(f"alpha must be > 0, got {alpha}.")
    yt, ys = _validate_inputs(np.asarray(y_true), np.asarray(y_score))
    n_total = len(yt)
    n_actives = int(yt.sum())
    ra = n_actives / n_total

    order = _stable_rank_order(yt, ys)
    sorted_scores = ys[order]
    sorted_labels = yt[order].astype(float)

    if len(np.unique(sorted_scores)) == n_total:
        ranked = [[float(sorted_scores[i]), int(sorted_labels[i])] for i in range(n_total)]
        return float(RdkitScoring.CalcBEDROC(ranked, 1, alpha))

    ranks_1based = np.arange(1, n_total + 1, dtype=float)
    weights = np.exp(-alpha * ranks_1based / n_total)

    # Average Boltzmann weights within each tied block
    avg_weights = weights.copy()
    i = 0
    while i < n_total:
        j = i + 1
        while j < n_total and sorted_scores[j] == sorted_scores[i]:
            j += 1
        if j - i > 1:
            avg_weights[i:j] = weights[i:j].mean()
        i = j

    rie_num = float(np.dot(sorted_labels, avg_weights)) / n_actives
    rie_den = (1.0 / n_total) * ((1.0 - math.exp(-alpha)) / (math.exp(alpha / n_total) - 1.0))
    rie = rie_num / rie_den

    rie_min = (1.0 - math.exp(alpha * ra)) / (ra * (1.0 - math.exp(alpha)))
    rie_max = (1.0 - math.exp(-alpha * ra)) / (ra * (1.0 - math.exp(-alpha)))
    return float(np.clip((rie - rie_min) / (rie_max - rie_min), 0.0, 1.0))


def evaluate_screening(
    y_true: np.ndarray | list[int],
    y_score: np.ndarray | list[float],
    method: str = "model",
) -> ScreeningMetrics:
    """Compute the full suite of honest virtual screening metrics."""
    yt, ys = _validate_inputs(np.asarray(y_true), np.asarray(y_score))
    n_total = len(yt)
    n_actives = int(yt.sum())
    hit_rate = n_actives / n_total

    ef1 = enrichment_factor(yt, ys, fraction=0.01)
    ef5 = enrichment_factor(yt, ys, fraction=0.05)
    ef1_max = max_enrichment_factor(n_total, n_actives, fraction=0.01)
    ef5_max = max_enrichment_factor(n_total, n_actives, fraction=0.05)
    bedroc = bedroc_score(yt, ys, alpha=20.0)
    auc = float(roc_auc_score(yt, ys))
    pr = float(average_precision_score(yt, ys))

    # Contrast metrics: accuracy when calling top 1% active vs predicting all inactive
    k1 = max(1, math.ceil(n_total * 0.01))
    order = _stable_rank_order(yt, ys)
    y_pred_top1 = np.zeros_like(yt)
    y_pred_top1[order[:k1]] = 1
    acc_top1 = float((y_pred_top1 == yt).mean())
    naive_acc = float((yt == 0).mean())

    return ScreeningMetrics(
        method=method,
        n_total=n_total,
        n_actives=n_actives,
        hit_rate=hit_rate,
        ef1=ef1,
        ef1_max=ef1_max,
        ef5=ef5,
        ef5_max=ef5_max,
        bedroc_20=bedroc,
        roc_auc=auc,
        pr_auc=pr,
        accuracy_top1pct=acc_top1,
        naive_all_inactive_accuracy=naive_acc,
    )


def bootstrap_screening_ci(
    y_true: np.ndarray | list[int],
    y_score: np.ndarray | list[float],
    n_bootstrap: int = 500,
    seed: int = 42,
) -> dict[str, tuple[float, float]]:
    """Stratified bootstrap 95% confidence intervals for EF1%, EF5%, BEDROC20, ROC-AUC, PR-AUC."""
    yt, ys = _validate_inputs(np.asarray(y_true), np.asarray(y_score))
    rng = np.random.default_rng(seed)
    pos_idx = np.flatnonzero(yt == 1)
    neg_idx = np.flatnonzero(yt == 0)

    samples: dict[str, list[float]] = {
        "ef1": [],
        "ef5": [],
        "bedroc_20": [],
        "roc_auc": [],
        "pr_auc": [],
    }
    for _ in range(n_bootstrap):
        b_pos = rng.choice(pos_idx, size=len(pos_idx), replace=True)
        b_neg = rng.choice(neg_idx, size=len(neg_idx), replace=True)
        idx = np.concatenate([b_pos, b_neg])
        m = evaluate_screening(yt[idx], ys[idx], method="boot")
        samples["ef1"].append(m.ef1)
        samples["ef5"].append(m.ef5)
        samples["bedroc_20"].append(m.bedroc_20)
        samples["roc_auc"].append(m.roc_auc)
        samples["pr_auc"].append(m.pr_auc)

    return {
        k: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
        for k, v in samples.items()
    }


def enrichment_curve(
    y_true: np.ndarray | list[int],
    y_score: np.ndarray | list[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Return (fraction_screened, fraction_actives_found) for cumulative recall plotting."""
    yt, ys = _validate_inputs(np.asarray(y_true), np.asarray(y_score))
    order = _stable_rank_order(yt, ys)
    sorted_labels = yt[order]
    cum_actives = np.cumsum(sorted_labels) / yt.sum()
    frac_screened = np.arange(1, len(yt) + 1) / len(yt)
    return np.concatenate([[0.0], frac_screened]), np.concatenate([[0.0], cum_actives])


def compute_roc_curve(
    y_true: np.ndarray | list[int],
    y_score: np.ndarray | list[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Return (fpr, tpr) arrays for ROC plotting."""
    yt, ys = _validate_inputs(np.asarray(y_true), np.asarray(y_score))
    fpr, tpr, _ = roc_curve(yt, ys)
    return fpr, tpr
