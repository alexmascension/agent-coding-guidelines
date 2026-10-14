import numpy as np
import pytest

from recipes.scrnaseq._clustering import best_match_jaccard
from recipes.scrnaseq.clustering_stability import choose_resolution, clustering_stability, default_n_iter
from recipes.scrnaseq.pc_robustness import (
    elbow_index,
    gavish_donoho_rank,
    pc_robustness,
    propose_n_pcs,
    robustness_verdict,
)

N_SIGNAL_DIMS = 6
N_CLUSTERS = 4


@pytest.fixture(scope="module")
def blobs():
    """Four well-separated clusters in 6 signal dimensions plus 24 noise dimensions, PCA-like ordering."""
    rng = np.random.default_rng(1)
    centres = rng.normal(0, 6, size=(N_CLUSTERS, N_SIGNAL_DIMS))
    labels = np.repeat(np.arange(N_CLUSTERS), 150)
    signal = centres[labels] + rng.normal(0, 1, size=(labels.size, N_SIGNAL_DIMS))
    noise = rng.normal(0, 1, size=(labels.size, 24))
    data = np.hstack([signal, noise])
    data -= data.mean(axis=0)
    u, s, vt = np.linalg.svd(data, full_matrices=False)
    return {
        "embedding": u * s,
        "singular_values": s,
        "shape": data.shape,
        "labels": labels,
        "variance_ratio": s**2 / np.sum(s**2),
    }


def test_best_match_jaccard_is_one_for_relabelled_partition():
    reference = np.array([0, 0, 1, 1, 2, 2])
    relabelled = np.array([5, 5, 3, 3, 9, 9])
    assert best_match_jaccard(reference, relabelled) == {0: 1.0, 1: 1.0, 2: 1.0}


def test_best_match_jaccard_rejects_different_lengths():
    with pytest.raises(ValueError, match="differ in length"):
        best_match_jaccard(np.array([0, 1]), np.array([0]))


def test_elbow_detects_signal_dimensions(blobs):
    # 4 clusters span at most 3 informative directions; the elbow must fall within the signal block
    assert 2 <= elbow_index(blobs["variance_ratio"]) <= N_SIGNAL_DIMS


def test_gavish_donoho_rank_is_within_signal_dimensions(blobs):
    rank = gavish_donoho_rank(blobs["singular_values"], *blobs["shape"])
    assert 2 <= rank <= N_SIGNAL_DIMS


def test_gavish_donoho_requires_full_spectrum(blobs):
    with pytest.raises(ValueError, match="full spectrum"):
        gavish_donoho_rank(blobs["singular_values"][:10], *blobs["shape"])


def test_propose_n_pcs_takes_larger_criterion(blobs):
    proposal = propose_n_pcs(blobs["variance_ratio"], blobs["singular_values"], blobs["shape"])
    assert proposal["proposal"] == max(proposal["elbow"], proposal["gavish_donoho"])


def test_pc_robustness_is_robust_for_separated_clusters(blobs):
    per_n, per_cluster = pc_robustness(blobs["embedding"], n_pcs=8, resolution=0.3, delta=2, seed=0)
    verdict = robustness_verdict(per_n, per_cluster)
    assert per_n.height == 5
    assert verdict["robust"], verdict


def test_pc_robustness_rejects_too_few_columns(blobs):
    with pytest.raises(ValueError, match="need"):
        pc_robustness(blobs["embedding"][:, :8], n_pcs=8, resolution=0.3, delta=2)


def test_default_n_iter_by_size():
    assert [default_n_iter(n) for n in (10_000, 100_000, 1_000_000)] == [10, 5, 1]


def test_clustering_stability_marks_true_clusters_stable(blobs):
    table = clustering_stability(blobs["embedding"], n_dims=8, resolutions=[0.3], n_iter=3, seed=0)
    assert set(table["status"].to_list()) == {"stable"}
    assert choose_resolution(table) == 0.3


def test_choose_resolution_returns_none_when_nothing_is_stable(blobs):
    table = clustering_stability(blobs["embedding"], n_dims=8, resolutions=[0.3], n_iter=2, seed=0)
    table = table.with_columns(status=np.array(["doubtful"] * table.height))
    assert choose_resolution(table) is None


def test_pc_robustness_flags_structureless_data():
    # pure noise has no reproducible partition: clusterings must change with the number of PCs
    rng = np.random.default_rng(2)
    noise = rng.normal(size=(500, 30))
    per_n, per_cluster = pc_robustness(noise, n_pcs=10, resolution=1.0, delta=3, seed=0)
    assert not robustness_verdict(per_n, per_cluster)["robust"]
