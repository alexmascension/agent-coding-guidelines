"""Recipe: cluster annotation from a marker dictionary with decoupler ULM, aggregated per cluster.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §9 (annotation and marker genes), first (automatic) step.

Procedure
---------
1. ``marker_network``: build the decoupler network from resolved marker IDs
   (``recipes.common.gene_ids``); positive
   markers weight +1, optional negative markers weight -1.
2. ``score_cells``: ULM (univariate linear model) score per cell and cell type on log-normalized
   expression.
3. ``annotate_clusters``: per cell, the label with the highest score; per cluster, the majority
   label, the fraction of
   cells that agree with it, and the median score of that label. A cluster is ``confident`` if
   agreement is at least
   ``MIN_AGREEMENT`` and the median score is positive; otherwise it is labelled ``unassigned`` for
   manual review.
   Cells whose top score is shared by several types are labelled ``tie``; cells with no score,
   ``unscored``.
   ``confident`` is a heuristic flag for the review, not a calibrated probability.

The result is a proposal for the manual review and verification steps, not a final annotation.
"""

import logging

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
import polars as pl

logger = logging.getLogger(__name__)

MIN_AGREEMENT = 0.5  # project default: majority of the cluster's cells share the label
MIN_MARKERS = 3  # decoupler ``tmin``: cell types with fewer markers present are not scored
UNASSIGNED = "unassigned"
TIE = "tie"  # per-cell label when the top score is shared by several cell types
UNSCORED = "unscored"  # per-cell label when no cell type could be scored


def marker_network(
    positive: dict[str, list[str]], negative: dict[str, list[str]] | None = None
) -> pd.DataFrame:
    """Return a decoupler network (``source``, ``target``, ``weight``) from marker ID dictionaries.

    Raises
    ------
    ValueError
        If a marker is both positive and negative for the same cell type.
    """
    for cell_type, genes in (negative or {}).items():
        both = sorted(set(genes) & set(positive.get(cell_type, [])))
        if both:
            raise ValueError(
                f"{cell_type}: markers listed as both positive and negative: {both[:5]}"
            )
    rows = [(cell_type, gene, 1.0) for cell_type, genes in positive.items() for gene in genes]
    for cell_type, genes in (negative or {}).items():
        rows += [(cell_type, gene, -1.0) for gene in genes]
    return pd.DataFrame(rows, columns=["source", "target", "weight"]).drop_duplicates(
        ["source", "target"]
    )


def score_cells(
    adata: ad.AnnData, net: pd.DataFrame, layer: str = "log1p_norm", tmin: int = MIN_MARKERS
) -> pd.DataFrame:
    """ULM scores per cell (rows) and cell type (columns) on ``layer``.

    Returns
    -------
    pd.DataFrame
        Cells x cell types; cell types with fewer than ``tmin`` markers in ``var_names`` are absent
        (logged).
    """
    tmp = ad.AnnData(
        X=adata.layers[layer],
        obs=pd.DataFrame(index=adata.obs_names),
        var=pd.DataFrame(index=adata.var_names),
    )
    dc.mt.ulm(data=tmp, net=net, tmin=tmin)
    scores = tmp.obsm["score_ulm"]
    scores = (
        scores if isinstance(scores, pd.DataFrame) else pd.DataFrame(scores, index=adata.obs_names)
    )
    dropped = sorted(set(net["source"]) - set(scores.columns))
    if dropped:
        logger.warning(f"cell types not scored (fewer than {tmin} markers present): {dropped}")
    return scores


def annotate_clusters(
    scores: pd.DataFrame, clusters: pd.Series, min_agreement: float = MIN_AGREEMENT
) -> pl.DataFrame:
    """Majority label per cluster from per-cell ULM scores.

    Parameters
    ----------
    scores : pd.DataFrame
        Output of ``score_cells``.
    clusters : pd.Series
        Cluster label per cell, same index as ``scores``.
    min_agreement : float
        Minimum fraction of cells whose top label equals the cluster label.

    Returns
    -------
    pl.DataFrame
        Per cluster: ``cluster``, ``n_cells``, ``label`` (or ``unassigned``), ``majority_label``,
        ``agreement``,
        ``median_score``, ``second_label``, ``second_fraction``, ``confident``.
    """
    if not clusters.index.equals(scores.index):
        raise ValueError("clusters and scores must have the same cell index, in the same order")
    values = scores.to_numpy(dtype=float)
    unscored = ~np.isfinite(values).any(axis=1)
    filled = np.where(np.isfinite(values), values, -np.inf)
    best = filled.max(axis=1, keepdims=True)
    tie = (filled == best).sum(axis=1) > 1
    top = pd.Series(
        np.asarray(scores.columns)[filled.argmax(axis=1)], index=scores.index, dtype=object
    )
    top[tie] = TIE
    top[unscored] = UNSCORED
    rows = []
    for cluster, idx in clusters.astype(str).groupby(clusters.astype(str)).groups.items():
        counts = top.loc[idx].value_counts(normalize=True).sort_index(kind="stable")
        counts = counts.sort_values(
            ascending=False, kind="stable"
        )  # ties between labels broken alphabetically
        majority, agreement = counts.index[0], float(counts.iloc[0])
        if majority in (TIE, UNSCORED):
            median_score, confident = float("nan"), False
        else:
            median_score = float(np.nanmedian(scores.loc[idx, majority]))
            confident = agreement >= min_agreement and median_score > 0
        rows.append(
            {
                "cluster": cluster,
                "n_cells": len(idx),
                "label": majority if confident else UNASSIGNED,
                "majority_label": majority,
                "agreement": agreement,
                "median_score": median_score,
                "second_label": counts.index[1] if counts.size > 1 else None,
                "second_fraction": float(counts.iloc[1]) if counts.size > 1 else 0.0,
                "confident": confident,
            }
        )
    return pl.DataFrame(rows).sort("cluster")
