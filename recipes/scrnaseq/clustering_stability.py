"""Recipe: choose the Leiden resolution by subsampling stability.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §8 (clustering).

Procedure
---------
For each resolution in a grid:
1. Cluster all cells (reference clustering).
2. ``n_iter`` times: subsample ``frac`` of the cells, rebuild the neighbour graph on the subsample from the same
   embedding and parameters (different seed per iteration), re-cluster.
3. For each reference cluster, Jaccard index between its subsampled cells and the best-matching subsample cluster.
4. Summarise per cluster: median and minimum Jaccard over iterations, and status (Hennig 2007 reading):
   median < 0.5 ``dissolved``; < 0.75 ``doubtful``; otherwise ``stable``.
``choose_resolution`` returns the highest resolution at which all (selected) clusters are stable.

Defaults
--------
* ``n_iter`` from ``default_n_iter``: 10 up to 50,000 cells, 5 up to 200,000, 1 above (indicative only).
* ``frac = 0.8``.

Typical use
-----------
>>> table = clustering_stability(adata.obsm["X_pca"], n_dims=30, resolutions=[0.2, 0.4, 0.6, 0.8, 1.0])
>>> best = choose_resolution(table)
"""

import logging

import numpy as np
import polars as pl

from recipes.scrnaseq._clustering import best_match_jaccard, leiden_labels

logger = logging.getLogger(__name__)

SUBSAMPLE_FRACTION = 0.8
JACCARD_DISSOLVED = 0.5  # Hennig 2007: clusters below this are not reproduced
JACCARD_STABLE = 0.75  # Hennig 2007: clusters at or above this are stable
N_ITER_BY_SIZE = ((50_000, 10), (200_000, 5))  # (max cells, iterations); above the last bound: 1 iteration


def default_n_iter(n_cells: int) -> int:
    """Return the number of subsampling iterations for a dataset size (project default, see module docstring)."""
    for max_cells, n_iter in N_ITER_BY_SIZE:
        if n_cells <= max_cells:
            return n_iter
    return 1


def _status(median_jaccard: float) -> str:
    if median_jaccard < JACCARD_DISSOLVED:
        return "dissolved"
    if median_jaccard < JACCARD_STABLE:
        return "doubtful"
    return "stable"


def clustering_stability(
    embedding: np.ndarray,
    n_dims: int,
    resolutions: list[float],
    n_neighbors: int = 15,
    frac: float = SUBSAMPLE_FRACTION,
    n_iter: int | None = None,
    seed: int = 0,
) -> pl.DataFrame:
    """Subsampling stability of Leiden clusters for each resolution.

    Parameters
    ----------
    embedding : np.ndarray
        Cells x components used for clustering (PCA or integrated embedding).
    n_dims : int
        Number of leading components used for the neighbour graph (the value chosen in §6).
    resolutions : list[float]
        Resolutions to evaluate.
    n_neighbors : int
        Neighbours of the k-NN graph, as in the main analysis.
    frac : float
        Fraction of cells kept in each subsample.
    n_iter : int | None
        Iterations per resolution; ``None`` uses ``default_n_iter``.
    seed : int
        Master seed; iteration ``i`` uses ``seed + 1 + i`` for subsampling and clustering.

    Returns
    -------
    pl.DataFrame
        One row per resolution x reference cluster: ``resolution`` (Float64), ``cluster`` (Int64), ``n_cells``
        (Int64), ``median_jaccard`` and ``min_jaccard`` (Float64), ``n_iter`` (Int64), ``status`` (Utf8).

    Raises
    ------
    ValueError
        If ``frac`` is not in (0, 1).
    """
    if not 0 < frac < 1:
        raise ValueError(f"frac must be in (0, 1), got {frac}")
    n_cells = embedding.shape[0]
    n_iter = default_n_iter(n_cells) if n_iter is None else n_iter
    if n_iter == 1:
        logger.warning(f"n_iter=1 for {n_cells} cells: stability is indicative only")
    n_sub = int(round(frac * n_cells))
    rows = []
    for resolution in resolutions:
        reference = leiden_labels(embedding, n_dims, resolution, n_neighbors, seed)
        jaccards: dict[int, list[float]] = {int(c): [] for c in np.unique(reference)}
        for i in range(n_iter):
            rng = np.random.default_rng(seed + 1 + i)
            idx = np.sort(rng.choice(n_cells, size=n_sub, replace=False))
            sub_labels = leiden_labels(embedding[idx], n_dims, resolution, n_neighbors, seed + 1 + i)
            for cluster, jac in best_match_jaccard(reference[idx], sub_labels).items():
                jaccards[cluster].append(jac)
        sizes = dict(zip(*np.unique(reference, return_counts=True), strict=True))
        for cluster, values in jaccards.items():
            # a reference cluster absent from every subsample cannot be assessed; it counts as dissolved
            median = float(np.median(values)) if values else 0.0
            rows.append(
                {
                    "resolution": float(resolution),
                    "cluster": cluster,
                    "n_cells": int(sizes[cluster]),
                    "median_jaccard": median,
                    "min_jaccard": float(min(values)) if values else 0.0,
                    "n_iter": n_iter,
                    "status": _status(median),
                }
            )
        logger.info(f"resolution {resolution:.2f}: {len(jaccards)} clusters")
    return pl.DataFrame(rows)


def choose_resolution(
    table: pl.DataFrame,
    clusters_of_interest: dict[float, list[int]] | None = None,
) -> float | None:
    """Highest resolution at which every (selected) cluster is ``stable``.

    Parameters
    ----------
    table : pl.DataFrame
        Output of ``clustering_stability``.
    clusters_of_interest : dict[float, list[int]] | None
        Optional, per resolution, the clusters that must be stable (labels differ between resolutions).

    Returns
    -------
    float | None
        Chosen resolution, or None if no resolution has all selected clusters stable.
    """
    candidates = []
    for (resolution,), group in table.group_by("resolution", maintain_order=True):
        if clusters_of_interest is not None and resolution in clusters_of_interest:
            group = group.filter(pl.col("cluster").is_in(clusters_of_interest[resolution]))
        if group.height and (group["status"] == "stable").all():
            candidates.append(resolution)
    return max(candidates) if candidates else None
