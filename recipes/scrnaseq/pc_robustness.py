"""Recipe: choose the number of principal components and check clustering robustness to it.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §6 (principal components).

Procedure
---------
1. ``propose_n_pcs``: proposal from the elbow of the explained-variance curve and, when the full
   singular-value
   spectrum is available, the Gavish-Donoho hard threshold; the larger of the two is proposed.
2. ``pc_robustness``: Leiden clustering at a fixed resolution for every number of PCs in ``p -
   delta .. p + delta``,
   compared with the clustering at ``p`` (ARI overall, best-match Jaccard per reference cluster).
3. ``robustness_verdict``: ``p`` is robust if ARI >= ``MIN_ARI`` for every tested value and every
   reference cluster
   keeps Jaccard >= ``MIN_JACCARD``.

Adapting
--------
* Pass the embedding actually used for clustering (``X_pca``, or an integrated one such as
  ``X_pca_harmony``).
* Keep ``n_neighbors``, ``resolution`` and ``seed`` identical to the main analysis.
* For large datasets use ``step=2``.

Typical use
-----------
>>> proposal = propose_n_pcs(adata.uns["pca"]["variance_ratio"])
>>> per_n, per_cluster = pc_robustness(adata.obsm["X_pca"], proposal["proposal"], resolution=0.8)
>>> verdict = robustness_verdict(per_n, per_cluster)
"""

import logging

import numpy as np
import polars as pl
from sklearn.metrics import adjusted_rand_score

from recipes.scrnaseq._clustering import best_match_jaccard, leiden_labels

logger = logging.getLogger(__name__)

MIN_ARI = (
    0.9  # project default: clusterings at neighbouring numbers of PCs should be nearly identical
)
MIN_JACCARD = (
    0.75  # stability reading of Hennig 2007 (clusterboot), as in SCRNASEQ_GUIDELINES.md §8
)
MAX_DENSE_ENTRIES = (
    2e8  # guard for the full SVD needed by Gavish-Donoho (cells x genes, dense float64)
)


def elbow_index(variance_ratio: np.ndarray) -> int:
    """Return the number of components at the elbow: maximum distance to the first-last chord.

    Parameters
    ----------
    variance_ratio : np.ndarray
        Explained variance ratio per component, in decreasing order.

    Returns
    -------
    int
        Number of components (1-based) at the elbow.
    """
    y = np.asarray(variance_ratio, dtype=float)
    if y.size < 3:
        raise ValueError(f"need at least 3 components, got {y.size}")
    x = np.arange(y.size, dtype=float)
    # normalise both axes so that the distance does not depend on units
    xn = (x - x[0]) / (x[-1] - x[0])
    yn = (y - y[-1]) / (y[0] - y[-1])
    # chord from (0, 1) to (1, 0): distance proportional to |xn + yn - 1|
    distance = np.abs(xn + yn - 1.0)
    return int(np.argmax(distance)) + 1


def gavish_donoho_rank(singular_values: np.ndarray, n_rows: int, n_cols: int) -> int:
    """Count singular values above the Gavish-Donoho optimal hard threshold (unknown noise level).

    Uses the approximation omega(beta) = 0.56 beta^3 - 0.95 beta^2 + 1.82 beta + 1.43 and the
    threshold
    omega(beta) * median(singular values) [Gavish 2014]. Requires the FULL spectrum (min(n_rows,
    n_cols) values).

    Parameters
    ----------
    singular_values : np.ndarray
        All singular values of the centred (and scaled, if used for PCA) matrix.
    n_rows, n_cols : int
        Matrix dimensions.

    Returns
    -------
    int
        Estimated number of signal components.
    """
    sv = np.sort(np.asarray(singular_values, dtype=float))[::-1]
    expected = min(n_rows, n_cols)
    if sv.size != expected:
        raise ValueError(f"need the full spectrum ({expected} values), got {sv.size}")
    beta = min(n_rows, n_cols) / max(n_rows, n_cols)
    omega = 0.56 * beta**3 - 0.95 * beta**2 + 1.82 * beta + 1.43
    return int(np.count_nonzero(sv > omega * np.median(sv)))


def full_singular_values(matrix: np.ndarray) -> np.ndarray:
    """All singular values of a dense, already centred/scaled matrix (cells x selected genes).

    Raises
    ------
    MemoryError
        If the matrix exceeds ``MAX_DENSE_ENTRIES``; then rely on the elbow only and say so.
    """
    if matrix.shape[0] * matrix.shape[1] > MAX_DENSE_ENTRIES:
        raise MemoryError(
            f"matrix {matrix.shape} too large for a full SVD; use the elbow criterion only"
        )
    return np.linalg.svd(np.asarray(matrix, dtype=float), compute_uv=False)


def propose_n_pcs(
    variance_ratio: np.ndarray,
    singular_values: np.ndarray | None = None,
    shape: tuple[int, int] | None = None,
    delta: int = 5,
) -> dict[str, int | str | None]:
    """Propose a number of PCs from the elbow and, if available, the Gavish-Donoho threshold.

    Parameters
    ----------
    variance_ratio : np.ndarray
        ``adata.uns["pca"]["variance_ratio"]``.
    singular_values : np.ndarray | None
        Full singular-value spectrum (see ``full_singular_values``); optional.
    shape : tuple[int, int] | None
        Shape of the matrix the spectrum comes from; required with ``singular_values``.
    delta : int
        Half-width of the robustness range that will be tested (``pc_robustness``).

    Returns
    -------
    dict[str, int | str | None]
        ``elbow``, ``gavish_donoho`` (None if not computed), ``proposal`` (the larger of the two,
        never truncated to
        the computed PCs), ``n_available``, ``n_required`` (``proposal + delta``) and ``status``:
        ``ok`` or
        ``recompute_pca`` when more components must be computed.
    """
    elbow = elbow_index(variance_ratio)
    gd = None
    if singular_values is not None:
        if shape is None:
            raise ValueError("shape is required when singular_values are given")
        gd = gavish_donoho_rank(singular_values, *shape)
    proposal = max(elbow, gd) if gd is not None else elbow
    n_available = len(variance_ratio)
    status = "ok" if proposal + delta <= n_available else "recompute_pca"
    if status != "ok":
        logger.warning(
            f"proposal {proposal} + delta {delta} exceeds the {n_available} computed PCs; "
            f"recompute PCA with at least {proposal + delta} components for the robustness test"
        )
    logger.info(f"PC proposal: elbow={elbow}, gavish_donoho={gd} -> {proposal}")
    return {
        "elbow": elbow,
        "gavish_donoho": gd,
        "proposal": proposal,
        "n_available": n_available,
        "n_required": proposal + delta,
        "status": status,
    }


def pc_robustness(
    embedding: np.ndarray,
    n_pcs: int,
    resolution: float,
    n_neighbors: int = 15,
    delta: int = 5,
    step: int = 1,
    seed: int = 0,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Compare clusterings at ``n_pcs - delta .. n_pcs + delta`` with the clustering at ``n_pcs``.

    Parameters
    ----------
    embedding : np.ndarray
        Cells x components; must have at least ``n_pcs + delta`` columns.
    n_pcs : int
        Proposed number of PCs (reference).
    resolution : float
        Leiden resolution of the main analysis.
    n_neighbors : int
        Neighbours of the k-NN graph of the main analysis.
    delta : int
        Half-width of the tested range.
    step : int
        Step between tested values (1 by default; 2 for large datasets).
    seed : int
        Seed shared by all clusterings, so that differences come from the number of PCs only.

    Returns
    -------
    per_n : pl.DataFrame
        One row per tested number of PCs: ``n_pcs`` (Int64), ``n_clusters`` (Int64), ``ari``
        (Float64).
    per_cluster : pl.DataFrame
        One row per tested number x reference cluster: ``n_pcs``, ``ref_cluster`` (Int64),
        ``ref_size`` (Int64),
        ``jaccard`` (Float64).

    Raises
    ------
    ValueError
        If the embedding has fewer than ``n_pcs + delta`` columns.
    """
    if embedding.shape[1] < n_pcs + delta:
        raise ValueError(
            f"embedding has {embedding.shape[1]} columns; need {n_pcs + delta}: "
            "recompute PCA with more components"
        )
    reference = leiden_labels(embedding, n_pcs, resolution, n_neighbors, seed)
    ref_sizes = {
        int(k): int(v) for k, v in zip(*np.unique(reference, return_counts=True), strict=True)
    }
    tested = sorted({n for n in range(n_pcs - delta, n_pcs + delta + 1, step) if n >= 2} | {n_pcs})
    rows_n, rows_c = [], []
    for n in tested:
        labels = (
            reference if n == n_pcs else leiden_labels(embedding, n, resolution, n_neighbors, seed)
        )
        rows_n.append(
            {
                "n_pcs": n,
                "n_clusters": int(np.unique(labels).size),
                "ari": float(adjusted_rand_score(reference, labels)),
            }
        )
        for cluster, jac in best_match_jaccard(reference, labels).items():
            rows_c.append(
                {
                    "n_pcs": n,
                    "ref_cluster": cluster,
                    "ref_size": ref_sizes[cluster],
                    "jaccard": jac,
                }
            )
    return pl.DataFrame(rows_n), pl.DataFrame(rows_c)


def robustness_verdict(
    per_n: pl.DataFrame,
    per_cluster: pl.DataFrame,
    clusters_of_interest: list[int] | None = None,
    min_ari: float = MIN_ARI,
    min_jaccard: float = MIN_JACCARD,
) -> dict[str, object]:
    """Summarise ``pc_robustness`` output into a decision.

    Returns
    -------
    dict[str, object]
        ``robust`` (bool), ``min_ari`` observed, and ``unstable_clusters``: reference clusters
        whose Jaccard fell
        below ``min_jaccard`` at any tested number of PCs.
    """
    selected = per_cluster
    if clusters_of_interest is not None:
        selected = per_cluster.filter(pl.col("ref_cluster").is_in(clusters_of_interest))
    unstable = (
        selected.group_by("ref_cluster", maintain_order=True)
        .agg(pl.col("jaccard").min().alias("min_jaccard"))
        .filter(pl.col("min_jaccard") < min_jaccard)
        .sort("ref_cluster")["ref_cluster"]
        .to_list()
    )
    observed_min_ari = float(per_n["ari"].min())
    return {
        "robust": observed_min_ari >= min_ari and not unstable,
        "min_ari": observed_min_ari,
        "unstable_clusters": unstable,
    }
