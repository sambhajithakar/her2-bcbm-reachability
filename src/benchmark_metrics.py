#!/usr/bin/env python3
"""Ranking metrics and inference for the leakage-controlled benchmark.

This module is complete and tested. It is not yet wired to any system output:
the agentic workflow has not been run, and the literature-volume counts have
not been collected. Running `python3 benchmark_metrics.py --selftest` checks
the implementations against hand-computed cases.

Primary endpoint (OSF decision 22): precision@20 of the agentic workflow
against the literature-volume baseline, compared by paired bootstrap.
"""
import argparse
import math

import numpy as np

RNG_DEFAULT = 20261013


def precision_at_k(ranked, relevant, k):
    if k <= 0:
        raise ValueError("k must be positive")
    top = ranked[:k]
    return sum(1 for item in top if item in relevant) / k


def recall_at_k(ranked, relevant, k):
    if not relevant:
        return float("nan")
    top = ranked[:k]
    return sum(1 for item in top if item in relevant) / len(relevant)


def dcg(gains):
    return sum(gain / math.log2(i + 2) for i, gain in enumerate(gains))


def ndcg_at_k(ranked, relevant, k):
    gains = [1.0 if item in relevant else 0.0 for item in ranked[:k]]
    ideal = [1.0] * min(len(relevant), k)
    denominator = dcg(ideal)
    return dcg(gains) / denominator if denominator else float("nan")


def mrr(ranked, relevant):
    for i, item in enumerate(ranked):
        if item in relevant:
            return 1.0 / (i + 1)
    return 0.0


def bootstrap_ci(ranked, relevant, metric, k=None, n_boot=2000,
                 seed=RNG_DEFAULT, alpha=0.05):
    """Percentile CI by resampling the ground-truth set, which is the quantity
    that is small and uncertain here. The ranking is fixed by the system."""
    rng = np.random.default_rng(seed)
    pool = list(relevant)
    if not pool:
        return float("nan"), float("nan")
    draws = []
    for _ in range(n_boot):
        sample = set(rng.choice(pool, size=len(pool), replace=True).tolist())
        draws.append(metric(ranked, sample, k) if k is not None
                     else metric(ranked, sample))
    lo, hi = np.quantile(draws, [alpha / 2, 1 - alpha / 2])
    return float(lo), float(hi)


def paired_bootstrap(ranked_a, ranked_b, relevant, metric, k=None,
                     n_boot=2000, seed=RNG_DEFAULT):
    """Two-sided paired bootstrap on the difference metric(A) - metric(B),
    resampling the ground-truth set. Returns (observed difference, p)."""
    rng = np.random.default_rng(seed)
    pool = list(relevant)
    call = (lambda r, s: metric(r, s, k)) if k is not None else metric
    observed = call(ranked_a, set(pool)) - call(ranked_b, set(pool))
    diffs = []
    for _ in range(n_boot):
        sample = set(rng.choice(pool, size=len(pool), replace=True).tolist())
        diffs.append(call(ranked_a, sample) - call(ranked_b, sample))
    diffs = np.asarray(diffs)
    centred = diffs - diffs.mean()
    p = float((np.abs(centred) >= abs(observed)).mean())
    return float(observed), p


def permutation_test(ranked_a, ranked_b, relevant, metric, k=None,
                     n_perm=5000, seed=RNG_DEFAULT):
    """Exchange the two rankings item-wise under the null that they are
    equivalent orderings of the same candidate pool."""
    rng = np.random.default_rng(seed)
    call = (lambda r, s: metric(r, s, k)) if k is not None else metric
    observed = call(ranked_a, relevant) - call(ranked_b, relevant)
    count = 0
    for _ in range(n_perm):
        swap = rng.random(len(ranked_a)) < 0.5
        perm_a = [b if s else a for a, b, s in zip(ranked_a, ranked_b, swap)]
        perm_b = [a if s else b for a, b, s in zip(ranked_a, ranked_b, swap)]
        if abs(call(perm_a, relevant) - call(perm_b, relevant)) >= abs(observed):
            count += 1
    return float(observed), (count + 1) / (n_perm + 1)


def brier_score(confidences, outcomes):
    """Calibration of the system's own stated confidence."""
    pairs = list(zip(confidences, outcomes))
    if not pairs:
        return float("nan")
    return sum((c - o) ** 2 for c, o in pairs) / len(pairs)


def expected_calibration_error(confidences, outcomes, n_bins=10):
    confidences = np.asarray(confidences, dtype=float)
    outcomes = np.asarray(outcomes, dtype=float)
    if confidences.size == 0:
        return float("nan")
    edges = np.linspace(0, 1, n_bins + 1)
    error = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (confidences > lo) & (confidences <= hi)
        if not mask.any():
            continue
        error += mask.mean() * abs(outcomes[mask].mean() - confidences[mask].mean())
    return float(error)


def _selftest():
    ranked = ["A", "B", "C", "D", "E"]
    relevant = {"B", "E", "Z"}

    assert precision_at_k(ranked, relevant, 2) == 0.5, "B is 1 of the top 2"
    assert precision_at_k(ranked, relevant, 5) == 0.4, "B and E are 2 of 5"
    assert abs(recall_at_k(ranked, relevant, 5) - 2 / 3) < 1e-12, "2 of 3 relevant found"
    assert mrr(ranked, relevant) == 0.5, "first hit at rank 2"
    assert mrr(["X", "Y"], relevant) == 0.0, "no hit"

    # nDCG by hand: gains 0,1,0,0,1 -> 1/log2(3) + 1/log2(6);
    # ideal has 3 relevant but k=5 caps at 3 -> 1/log2(2)+1/log2(3)+1/log2(4)
    got = ndcg_at_k(ranked, relevant, 5)
    want = (1 / math.log2(3) + 1 / math.log2(6)) / (
        1 / math.log2(2) + 1 / math.log2(3) + 1 / math.log2(4))
    assert abs(got - want) < 1e-12, f"nDCG {got} != {want}"

    perfect = ["B", "E", "Z", "A", "C"]
    assert abs(ndcg_at_k(perfect, relevant, 3) - 1.0) < 1e-12, "perfect ranking is 1.0"

    lo, hi = bootstrap_ci(ranked, relevant, precision_at_k, k=5)
    assert 0.0 <= lo <= hi <= 1.0, "CI must be ordered and in range"

    diff, p = paired_bootstrap(perfect, ranked, relevant, precision_at_k, k=5)
    assert diff > 0, "the perfect ranking must score higher"
    assert 0.0 <= p <= 1.0

    diff, p = paired_bootstrap(ranked, ranked, relevant, precision_at_k, k=5)
    assert diff == 0.0 and p == 1.0, "identical rankings differ by exactly nothing"

    assert abs(brier_score([1.0, 0.0], [1.0, 0.0])) < 1e-12, "perfect forecasts"
    assert abs(brier_score([0.5, 0.5], [1.0, 0.0]) - 0.25) < 1e-12

    assert abs(expected_calibration_error([0.9, 0.9], [1.0, 1.0])
               - 0.1) < 1e-12, "confident and right, off by 0.1"

    print("all self-tests passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        _selftest()
    else:
        parser.print_help()
