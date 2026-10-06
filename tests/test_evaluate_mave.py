"""Metric helpers for MAVE evaluation (ties by expectation)."""

from tessera.evaluate_mave import auroc, avg_ranks, expected_hits_at_k, spearman


def test_auroc_counts_ties_as_half():
    assert auroc([(2,), (1,)], [1, 0]) == 1.0
    assert auroc([(1,), (1,)], [1, 0]) == 0.5


def test_expected_hits_split_the_boundary_tie_group():
    keys = [(3,), (2,), (2,), (2,), (1,)]
    labels = [1, 1, 0, 0, 1]
    assert expected_hits_at_k(keys, labels, 2) == 1 + 1 / 3


def test_avg_ranks_and_spearman():
    assert avg_ranks([(3,), (1,), (3,)]) == [1.5, 3.0, 1.5]
    assert abs(spearman([1, 2, 3], [10, 20, 30]) - 1.0) < 1e-12
