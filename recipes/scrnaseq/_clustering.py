"""Shared helpers for clustering-based recipes (neighbour graph, Leiden, cluster matching).

Used by ``pc_robustness.py`` and ``clustering_stability.py``. Guideline:
``SCRNASEQ_GUIDELINES.md`` §6 (principal components) and §8 (clustering).
"""

import hashlib
import logging
from collections.abc import Sequence
from typing import Literal

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp

logger = logging.getLogger(__name__)

LEIDEN_N_ITERATIONS = (
    2  # igraph flavour; 2 iterations is the scanpy recommendation for this flavour
)
KNN_DECIMALS = 6  # rounding of the embedding before the tree search (SC-CLUST-004)
ZERO_DISTANCE = 1e-12  # stands for 0 between duplicated cells: scanpy drops stored zeros


def tree_knn_distances(
    embedding: np.ndarray, n_neighbors: int = 15, decimals: int = KNN_DECIMALS
) -> sp.csr_matrix:
    """Exact k-NN distances that do not depend on the thread count or the CPU (``SC-CLUST-004``).

    A kd-tree (no BLAS, one thread) on the embedding in float64 rounded to ``decimals``. The result
    follows the convention of ``scanpy.pp.neighbors``, whose ``n_neighbors`` counts the cell
    itself: ``n_neighbors - 1`` neighbours per cell, the cell excluded. Exact zeros (duplicated
    cells) become ``ZERO_DISTANCE``, because scanpy drops stored zeros and the row would lose
    neighbours. Use as ``sc.pp.neighbors(adata, n_neighbors=n_neighbors, distances=graph)``.

    Parameters
    ----------
    embedding : np.ndarray
        Cells x components (e.g. the first PCs).
    n_neighbors : int
        Number of neighbours as in ``scanpy.pp.neighbors`` (the cell included).
    decimals : int
        Decimals kept before the search; absorbs differences of the last bits between BLAS builds.

    Returns
    -------
    sp.csr_matrix
        Cells x cells distance matrix with ``n_neighbors - 1`` stored entries per row.

    Raises
    ------
    ValueError
        If ``n_neighbors`` < 2 or there are not more cells than ``n_neighbors - 1``.
    """
    from sklearn.neighbors import NearestNeighbors

    if n_neighbors < 2:
        raise ValueError(f"n_neighbors={n_neighbors} must be at least 2")
    if embedding.shape[0] <= n_neighbors - 1:
        raise ValueError(f"{embedding.shape[0]} cells are too few for n_neighbors={n_neighbors}")
    points = np.round(np.asarray(embedding, dtype=np.float64), decimals)
    search = NearestNeighbors(n_neighbors=n_neighbors - 1, algorithm="kd_tree", n_jobs=1)
    graph = sp.csr_matrix(search.fit(points).kneighbors_graph(mode="distance"))
    graph.data[graph.data == 0] = ZERO_DISTANCE
    return graph


def labels_fingerprint(obs_names: Sequence[str], labels: Sequence) -> str:
    """SHA-256 of the barcode-to-label pairs sorted by barcode (``SC-CLUST-005``).

    Store it with any decision keyed by cluster label and refuse the decisions when the fingerprint
    of the current clustering differs.
    """
    names = [str(n) for n in obs_names]
    values = [str(v) for v in labels]
    if len(names) != len(values):
        raise ValueError(f"{len(names)} names but {len(values)} labels")
    pairs = sorted(zip(names, values, strict=True))
    return hashlib.sha256("".join(f"{n}\t{v}\n" for n, v in pairs).encode()).hexdigest()


def leiden_labels(
    embedding: np.ndarray,
    n_dims: int,
    resolution: float,
    n_neighbors: int = 15,
    seed: int = 0,
    knn: Literal["approximate", "tree"] = "approximate",
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
    knn : {"approximate", "tree"}
        ``approximate``: scanpy's default search (exploration); ``tree``: ``tree_knn_distances``,
        the same partition for any thread count and CPU (``SC-CLUST-004``).

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
    if knn == "tree":
        graph = tree_knn_distances(tmp.obsm["emb"], n_neighbors)
        sc.pp.neighbors(tmp, n_neighbors=n_neighbors, distances=graph)
    elif knn == "approximate":
        sc.pp.neighbors(tmp, n_neighbors=n_neighbors, use_rep="emb", random_state=seed)
    else:
        raise ValueError(f"knn must be 'approximate' or 'tree', not {knn!r}")
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
