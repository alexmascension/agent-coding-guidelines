"""Shared helpers for clustering-based recipes (neighbour graph, Leiden, cluster matching).

Used by ``pc_robustness.py`` and ``clustering_stability.py``. Guideline:
``SCRNASEQ_GUIDELINES.md`` §6 (principal components) and §8 (clustering).
"""

import logging

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp

logger = logging.getLogger(__name__)

LEIDEN_N_ITERATIONS = (
    2  # igraph flavour; 2 iterations is the scanpy recommendation for this flavour
)


def leiden_labels(
    embedding: np.ndarray,
    n_dims: int,
    resolution: float,
    n_neighbors: int = 15,
    seed: int = 0,
) -> np.ndarray:
    """Cluster cells with Leiden on a k-NN graph of the first ``n_dims`` embedding columns.

    Parameters
    ----------
    embedding : np.ndarray
        Cells x components (e.g. ``adata.obsm["X_pca"]`` or an integrated embedding such as
        Harmony).
    n_dims : int
        Number of leading components used to build the neighbour graph.
    resolution : float
        Leiden resolution.
    n_neighbors : int
        Number of neighbours of the k-NN graph.
    seed : int
        Random seed for the neighbour search and Leiden.

    Returns
    -------
    np.ndarray
        Integer cluster label per cell.

    Raises
    ------
    ValueError
        If ``n_dims`` exceeds the number of available components.
    """
    if n_dims > embedding.shape[1]:
        raise ValueError(f"n_dims={n_dims} exceeds available components ({embedding.shape[1]})")
    tmp = ad.AnnData(obs={"i": np.arange(embedding.shape[0])})
    tmp.obsm["emb"] = np.ascontiguousarray(embedding[:, :n_dims])
    sc.pp.neighbors(tmp, n_neighbors=n_neighbors, use_rep="emb", random_state=seed)
    sc.tl.leiden(
        tmp,
        resolution=resolution,
        flavor="igraph",
        n_iterations=LEIDEN_N_ITERATIONS,
        directed=False,
        random_state=seed,
    )
    return tmp.obs["leiden"].astype(int).to_numpy()


def best_match_jaccard(reference: np.ndarray, other: np.ndarray) -> dict[int, float]:
    """Jaccard index of each reference cluster with its best-matching cluster in another labelling.

    Both labellings refer to the same cells, in the same order. For a reference cluster A and a
    cluster B of the other
    labelling, J = |A & B| / (|A| + |B| - |A & B|); the best match maximises J. Computed from the
    contingency table,
    so memory grows with the number of cluster pairs, not with cells x clusters.

    Parameters
    ----------
    reference : np.ndarray
        Labels of the reference clustering.
    other : np.ndarray
        Labels of the clustering to compare, same length as ``reference``.

    Returns
    -------
    dict[int, float]
        Reference cluster label -> best Jaccard index.

    Raises
    ------
    ValueError
        If the two label arrays differ in length.
    """
    reference, other = np.asarray(reference), np.asarray(other)
    if reference.shape != other.shape:
        raise ValueError(f"label arrays differ in length: {reference.shape} vs {other.shape}")
    ref_levels, ref_codes = np.unique(reference, return_inverse=True)
    other_levels, other_codes = np.unique(other, return_inverse=True)
    table = sp.coo_matrix(
        (np.ones(reference.size), (ref_codes, other_codes)),
        shape=(ref_levels.size, other_levels.size),
    ).toarray()
    ref_sizes = table.sum(axis=1, keepdims=True)
    other_sizes = table.sum(axis=0, keepdims=True)
    jaccard = table / (ref_sizes + other_sizes - table)
    return {
        int(level): float(best)
        for level, best in zip(ref_levels, jaccard.max(axis=1), strict=True)
    }
