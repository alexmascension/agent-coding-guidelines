"""Recipe: cluster annotation from a marker dictionary with decoupler ULM, aggregated per cluster.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §9 (annotation and marker genes), first (automatic) step.

Procedure
---------
1. ``marker_network``: build the decoupler network from resolved marker IDs (``recipes.common.gene_ids``); positive
   markers weight +1, optional negative markers weight -1.
2. ``score_cells``: ULM (univariate linear model) score per cell and cell type on log-normalized expression.
3. ``annotate_clusters``: per cell, the label with the highest score; per cluster, the majority label, the fraction of
   cells that agree with it, and the median score of that label. A cluster is ``confident`` if agreement is at least
   ``MIN_AGREEMENT`` and the median score is positive; otherwise it is labelled ``unassigned`` for manual review.

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


def marker_network(positive: dict[str, list[str]], negative: dict[str, list[str]] | None = None) -> pd.DataFrame:
    """Return a decoupler network (``source``, ``target``, ``weight``) from marker ID dictionaries."""
    rows = [(cell_type, gene, 1.0) for cell_type, genes in positive.items() for gene in genes]
    for cell_type, genes in (negative or {}).items():
        rows += [(cell_type, gene, -1.0) for gene in genes]
    net = pd.DataFrame(rows, columns=["source", "target", "weight"]).drop_duplicates(["source", "target"])
    return net


def score_cells(
    adata: ad.AnnData, net: pd.DataFrame, layer: str = "log1p_norm", tmin: int = MIN_MARKERS
) -> pd.DataFrame:
    """ULM scores per cell (rows) and cell type (columns) on ``layer``.

    Returns
    -------
    pd.DataFrame
        Cells x cell types; cell types with fewer than ``tmin`` markers in ``var_names`` are absent (logged).
    """
    tmp = ad.AnnData(
        X=adata.layers[layer], obs=pd.DataFrame(index=adata.obs_names), var=pd.DataFrame(index=adata.var_names)
    )
    dc.mt.ulm(data=tmp, net=net, tmin=tmin)
    scores = tmp.obsm["score_ulm"]
    scores = scores if isinstance(scores, pd.DataFrame) else pd.DataFrame(scores, index=adata.obs_names)
    dropped = sorted(set(net["source"]) - set(scores.columns))
    if dropped:
        logger.warning(f"cell types not scored (fewer than {tmin} markers present): {dropped}")
    return scores


def annotate_clusters(scores: pd.DataFrame, clusters: pd.Series, min_agreement: float = MIN_AGREEMENT) -> pl.DataFrame:
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
        Per cluster: ``cluster``, ``n_cells``, ``label`` (or ``unassigned``), ``majority_label``, ``agreement``,
        ``median_score``, ``second_label``, ``second_fraction``, ``confident``.
    """
    top = scores.idxmax(axis=1)
    rows = []
    for cluster, idx in clusters.astype(str).groupby(clusters.astype(str)).groups.items():
        counts = top.loc[idx].value_counts(normalize=True)
        majority, agreement = counts.index[0], float(counts.iloc[0])
        median_score = float(np.median(scores.loc[idx, majority]))
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
