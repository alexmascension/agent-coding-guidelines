"""Recipe: compare cell types at matched depth (cells per unit and UMIs per cell).

Guideline: ``SCRNASEQ_GUIDELINES.md`` §10 (``SC-DE-008``): the number of DE genes per cell type
depends on the donors with that cell type, the cells per donor and the UMIs per cell, so cell types
are compared at matched depth, or by effect sizes.

``match_depth`` reduces a target cell type to the depth of a reference cell type in the same units:

1. cells: in each unit, a random subset of the target cells, as many as the reference cell type has
   in that unit (fewer if the target has fewer); units without reference cells lose the target
   cells too;
2. UMIs: binomial thinning of every count of the kept target cells with probability
   ``f = mean UMIs per reference cell / mean UMIs per kept target cell``, which keeps the expected
   proportions of genes in each cell and lowers its total to ``f`` times on average. Thinning
   matches the shallower sequencing of the reference in expectation; it does not change cell
   composition, ambient RNA or capture efficiency, which also differ between cell types.

Matching serves the comparison between cell types; the main analysis uses all the data. Counts
MUST be integers (raw or CellBender counts); SoupX- or decontX-corrected layers are not, unless
they were rounded explicitly and the rounding recorded. The default layer is resolved with
``recipes/scrnaseq/layers.py`` (``SC-DATA-004``, step ``"pseudobulk"``); a non-integer layer
raises an error that names it.
"""

import logging

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp

from recipes.scrnaseq.layers import input_layer, require_integer_counts

logger = logging.getLogger(__name__)

MATCH_MODES = ("cells", "umis")


def _integer_csr(matrix) -> sp.csr_matrix:
    """Return ``matrix`` as CSR with int64 data, refusing negative or non-integer values."""
    csr = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
    data = csr.data
    if np.any(data < 0) or np.any(data != np.round(data)):
        raise ValueError("counts must be non-negative integers (raw or CellBender counts)")
    out = csr.copy()
    out.data = np.round(data).astype(np.int64)
    return out


def match_depth(
    adata: ad.AnnData,
    unit_key: str,
    celltype_key: str,
    target: str,
    reference: str,
    layer: str | None = None,
    match: tuple[str, ...] = MATCH_MODES,
    seed: int = 0,
) -> tuple[ad.AnnData, pl.DataFrame, dict[str, float]]:
    """Reduce the cells of ``target`` to the cells per unit and UMIs per cell of ``reference``.

    Parameters
    ----------
    adata : ad.AnnData
        Cells after QC, with integer counts in ``layer``.
    unit_key : str
        ``obs`` column with the aggregation unit (donor or sample).
    celltype_key : str
        ``obs`` column with the annotation.
    target : str
        Cell type to reduce.
    reference : str
        Cell type whose depth is matched (the shallower one).
    layer : str | None
        Integer counts layer; None resolves it with ``input_layer(adata, "pseudobulk")``.
    match : tuple[str, ...]
        What to match: ``"cells"``, ``"umis"`` or both.
    seed : int
        Seed of the cell subsampling and of the thinning.

    Returns
    -------
    matched : ad.AnnData
        The kept target cells, with the matched counts in ``layer`` and in ``X``; ``obs`` and
        ``var`` copied from the input.
    per_unit : pl.DataFrame
        ``unit``, ``n_target``, ``n_reference``, ``n_kept``.
    summary : dict[str, float]
        ``umi_factor`` (1 when no thinning was applied), ``mean_umis_target_before``,
        ``mean_umis_target_after``, ``mean_umis_reference``, ``n_cells_before``,
        ``n_cells_after``.

    Raises
    ------
    ValueError
        If ``match`` is empty or unknown, a cell type is absent, or counts are not integers.
    """
    if not match or set(match) - set(MATCH_MODES):
        raise ValueError(f"match must be a non-empty subset of {MATCH_MODES}, got {match}")
    cell_types = adata.obs[celltype_key].astype(str)
    for cell_type in (target, reference):
        if not (cell_types == cell_type).any():
            raise ValueError(f"cell type {cell_type!r} not in {celltype_key!r}")
    layer = input_layer(adata, "pseudobulk") if layer is None else layer
    require_integer_counts(adata, layer)
    rng = np.random.default_rng(seed)
    counts = _integer_csr(adata.layers[layer])
    units = adata.obs[unit_key].astype(str).to_numpy()
    is_target = (cell_types == target).to_numpy()
    is_reference = (cell_types == reference).to_numpy()
    n_target = pd.Series(units[is_target]).value_counts()
    n_reference = pd.Series(units[is_reference]).value_counts()

    target_index = np.flatnonzero(is_target)
    if "cells" in match:
        kept, rows = [], []
        for unit in sorted(n_target.index):
            in_unit = target_index[units[target_index] == unit]
            k = min(in_unit.size, int(n_reference.get(unit, 0)))
            kept.append(rng.choice(in_unit, size=k, replace=False))
            rows.append(
                {
                    "unit": unit,
                    "n_target": int(in_unit.size),
                    "n_reference": int(n_reference.get(unit, 0)),
                    "n_kept": int(k),
                }
            )
        kept_index = np.sort(np.concatenate(kept)) if kept else np.array([], dtype=int)
    else:
        kept_index = target_index
        rows = [
            {
                "unit": unit,
                "n_target": int(n),
                "n_reference": int(n_reference.get(unit, 0)),
                "n_kept": int(n),
            }
            for unit, n in sorted(n_target.items())
        ]

    sub = counts[kept_index]
    before = np.asarray(sub.sum(axis=1)).ravel()
    mean_reference = float(np.asarray(counts[is_reference].sum(axis=1)).mean())
    mean_before = float(before.mean()) if before.size else float("nan")
    factor = 1.0
    if "umis" in match and before.size:
        factor = mean_reference / mean_before
        if factor >= 1:
            logger.info(
                f"{target} is not deeper than {reference} (factor {factor:.3f}); no thinning"
            )
            factor = 1.0
        else:
            sub = sub.copy()
            sub.data = rng.binomial(sub.data, factor)
            sub.eliminate_zeros()
    after = np.asarray(sub.sum(axis=1)).ravel()
    matched = ad.AnnData(
        X=sub.astype(np.float32),
        obs=adata.obs.iloc[kept_index].copy(),
        var=adata.var.copy(),
    )
    matched.layers[layer] = sub.astype(np.float32)
    summary = {
        "umi_factor": factor,
        "mean_umis_target_before": mean_before,
        "mean_umis_target_after": float(after.mean()) if after.size else float("nan"),
        "mean_umis_reference": mean_reference,
        "n_cells_before": int(is_target.sum()),
        "n_cells_after": int(kept_index.size),
    }
    return matched, pl.DataFrame(rows), summary
