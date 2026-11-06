"""Recipe: pseudobulk aggregation per unit and cell type, condition DE per cell type, diagnostics.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §10 (differential expression between conditions); DE steps
follow ``RNASEQ_BULK_GUIDELINES.md`` §5-8 through ``recipes/bulk/deseq2.py``.

Procedure
---------
1. ``aggregate_counts``: sum counts per (aggregation unit, cell type) and record the number of
   cells. The aggregation unit is the experimental unit of the comparison: a sample, or a donor
   only when each donor contributes to one level of every compared factor. When a donor has
   samples in several conditions or time points, aggregate per sample and model the donor
   (``SCRNASEQ_GUIDELINES.md`` §2). With ``sample_keys``, the cells of one physical sample of one
   donor MUST fall in a single unit; splitting them creates pseudo-replicates and raises an error
   (``SC-DE-006``). Missing unit identifiers, and distinct identifiers that become the same string
   (``1`` and ``"1"``), raise an error before any conversion (``recipes/scrnaseq/units.py``). The
   default layer is resolved with ``recipes/scrnaseq/layers.py`` (``SC-DATA-004``:
   ``counts_ambient_corrected`` if it exists, else ``counts``); non-integer counts raise an error
   that names the layer, and are never rounded here.
2. Units with fewer than ``min_cells`` cells (default ``MIN_CELLS``) are not used. The report lists
   every (unit, cell type) combination with status ``kept``, ``below_min_cells`` or ``absent`` (no
   cells). A cell type absent from a unit is not turned into a pseudobulk sample of zeros.
   ``min_cells_diagnostic`` (on an aggregation with ``min_cells=1``) shows where samples with few
   cells drift away from the rest, to fix the threshold before the DE results are seen
   (``SC-DE-002``).
3. ``pseudobulk_de``: per cell type, design-based gene filter and pydeseq2 with the project design;
   cell types where a compared group has fewer than ``MIN_SAMPLES_PER_GROUP`` units are reported
   as ``not_testable``. Adjusted p-values are given per cell type (``p_adjusted``, pydeseq2 with
   independent filtering within the cell type) and across all cell types (``p_adjusted_global``,
   BH over the raw p-values of all tested genes). The two families differ, also because
   independent filtering is applied per cell type; the project states which one it reports
   (``STATISTICS_GUIDELINES.md`` §7, multiple testing). Each gene also gets the Wald interval of
   its log fold change (``lfc_ci_low``, ``lfc_ci_high``), its rank by p-value within the cell type
   and, when significant, the quartile of its adjusted p-value among the significant genes
   (``significance_quartile``, 1 = most significant; ``SC-DE-010``).
4. ``detectability_report``: for the genes the question is about, whether each was tested, changed,
   shown not to change by more than a margin (equivalence test), or is inconclusive, with the
   expression, interval and units behind it (``SC-DE-009``).
5. ``null_split_de``: global-null diagnostic per cell type, by random splits of the control units
   or by permuted condition labels, through ``recipes/bulk/fdr_calibration.py`` (``SC-DE-007``).
6. ``power_inputs``: per cell type and group, units, cells per unit and UMIs per cell, the
   quantities that decide power and that comparisons between cell types report (``SC-DE-008``,
   ``SC-REP-006``).

Covariates keep their dtype: numeric covariates (age, RIN, postmortem interval) stay numeric and
enter the design as slopes; categorical ones stay categorical. Missing values MUST be resolved
before calling (an error is raised).
"""

import logging
import re

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests

from recipes.bulk.deseq2 import RESULT_COLUMNS, filter_genes, run_deseq2
from recipes.bulk.fdr_calibration import fdr_calibration
from recipes.scrnaseq.layers import input_layer, require_integer_counts
from recipes.scrnaseq.units import unit_ids, unit_metadata

logger = logging.getLogger(__name__)

MIN_CELLS = 10  # guide: 10-20 cells per pseudobulk sample
MIN_SAMPLES_PER_GROUP = (
    2  # computational minimum for a dispersion estimate, not a sample-size justification
)
ALPHA = 0.05  # FDR level of the DE calls
CI_LEVEL = 0.95  # level of the Wald interval of each log fold change
N_NULL_SPLITS = 20  # random splits (or permutations) of the global-null diagnostic
NULL_SPLIT_COLUMN = "null_split"  # label column created by null_split_de
NULL_SPLIT_LEVELS = ("split_a", "split_b")
DIAGNOSTIC_PCS = 2  # principal components used by min_cells_diagnostic
DIAGNOSTIC_GENES = 2000  # most variable genes (log CPM) used by min_cells_diagnostic
DIAGNOSTIC_MIN_CPM = 1.0  # genes below this mean CPM are not used by min_cells_diagnostic
CELL_BINS = (1, 5, 10, 20, 30, 50, 100)  # lower edges of the cells-per-unit bins of the diagnostic
REPORT_SCHEMA = {"cell_type": pl.Utf8, "unit": pl.Utf8, "n_cells": pl.Int64, "status": pl.Utf8}
DE_EXTRA_COLUMNS = [
    "cell_type",
    "p_adjusted_global",
    "lfc_ci_low",
    "lfc_ci_high",
    "rank",
    "significance_quartile",
]


def check_no_pseudo_replicates(obs: pd.DataFrame, unit_key: str, sample_keys: list[str]) -> None:
    """Raise if the cells of one physical sample are split across several aggregation units.

    Parameters
    ----------
    obs : pd.DataFrame
        Cell metadata.
    unit_key : str
        Column with the aggregation unit.
    sample_keys : list[str]
        Columns that together identify one physical sample of one donor: ``["sample_id"]``, or
        ``["sample_id", "donor_id"]`` for captures that pool several donors.

    Raises
    ------
    ValueError
        If a sample maps to more than one unit (pseudo-replicates, ``SC-DE-006``), a key is
        absent, or a key has missing values or values that collide as strings (``unit_ids``).
    """
    absent = [k for k in [unit_key, *sample_keys] if k not in obs.columns]
    if absent:
        raise ValueError(f"columns not in obs: {absent}")
    keys = {k: unit_ids(obs, k) for k in sample_keys}
    n_units = unit_ids(obs, unit_key).groupby([keys[k] for k in sample_keys]).nunique()
    split = n_units[n_units > 1]
    if not split.empty:
        examples = [str(i) for i in split.index[:5]]
        raise ValueError(
            f"cells of {split.size} sample(s) of {sample_keys} are split across several values of "
            f"{unit_key!r} (e.g. {examples}): these are pseudo-replicates, not experimental units "
            "(SCRNASEQ_GUIDELINES.md SC-DE-006)"
        )


def aggregate_counts(
    adata: ad.AnnData,
    unit_key: str,
    celltype_key: str,
    covariates: list[str],
    layer: str | None = None,
    min_cells: int = MIN_CELLS,
    sample_keys: list[str] | None = None,
) -> tuple[dict[str, tuple[pd.DataFrame, pd.DataFrame]], pl.DataFrame]:
    """Sum counts per aggregation unit and cell type.

    Parameters
    ----------
    adata : ad.AnnData
        Cells after QC and doublet removal; integer counts in ``layer``.
    unit_key : str
        ``obs`` column with the experimental unit (sample, or donor when appropriate; see module
        docstring). Missing identifiers and values that collide as strings are rejected.
    celltype_key : str
        ``obs`` column with the final annotation.
    covariates : list[str]
        Unit-level ``obs`` columns carried into the metadata; dtypes are preserved.
    layer : str | None
        Counts layer; None resolves it with ``input_layer(adata, "pseudobulk")``
        (``counts_ambient_corrected`` if it exists, else ``counts``; ``SC-DATA-004``). The layer
        MUST hold integer counts (``require_integer_counts``).
    min_cells : int
        Minimum cells per pseudobulk sample.
    sample_keys : list[str] | None
        Columns identifying one physical sample of one donor (``check_no_pseudo_replicates``);
        None skips the check.

    Returns
    -------
    pseudobulk : dict[str, tuple[pd.DataFrame, pd.DataFrame]]
        Cell type -> (counts units x genes, metadata units x [covariates..., ``n_cells``]), kept
        units only.
    report : pl.DataFrame
        Every (cell type, unit) combination: ``cell_type``, ``unit``, ``n_cells``, ``status``.

    Raises
    ------
    ValueError
        If unit identifiers or covariates are invalid (``unit_metadata``), cells of one sample
        fall in several units (``check_no_pseudo_replicates``), or the layer does not hold
        integer counts.
    """
    if sample_keys:
        check_no_pseudo_replicates(adata.obs, unit_key, sample_keys)
    meta_all = unit_metadata(adata.obs, unit_key, covariates)
    layer = input_layer(adata, "pseudobulk") if layer is None else layer
    require_integer_counts(adata, layer)
    units = adata.obs[unit_key].astype(str).to_numpy()
    cell_types = adata.obs[celltype_key].astype(str).to_numpy()
    counts = adata.layers[layer]
    counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)
    all_units = meta_all.index.to_numpy()

    pseudobulk, rows = {}, []
    for cell_type in np.unique(cell_types):
        in_type = np.flatnonzero(cell_types == cell_type)
        unit_index = pd.Index(all_units).get_indexer(units[in_type])
        indicator = sp.csr_matrix(
            (np.ones(in_type.size), (unit_index, np.arange(in_type.size))),
            shape=(all_units.size, in_type.size),
        )
        summed = (indicator @ counts[in_type]).toarray()
        n_cells = np.bincount(unit_index, minlength=all_units.size)
        kept = n_cells >= min_cells
        for unit, n, k in zip(all_units, n_cells, kept, strict=True):
            status = "kept" if k else ("absent" if n == 0 else "below_min_cells")
            rows.append(
                {"cell_type": cell_type, "unit": unit, "n_cells": int(n), "status": status}
            )
        if not kept.any():
            continue
        frame = pd.DataFrame(summed[kept], index=all_units[kept], columns=adata.var_names)
        meta = meta_all.loc[all_units[kept]].copy()
        meta["n_cells"] = n_cells[kept]
        pseudobulk[cell_type] = (frame, meta)
    return pseudobulk, pl.DataFrame(rows, schema=REPORT_SCHEMA)


def _rank_and_quartile(table: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """Add ``rank`` (by p-value) and ``significance_quartile`` within one cell type."""
    table = table.copy()
    table["rank"] = table["p_value"].rank(method="min", na_option="keep")
    significant = table["p_adjusted"] < alpha
    quartile = pd.Series(np.nan, index=table.index)
    n_significant = int(significant.sum())
    if n_significant:
        order = table.loc[significant, "p_adjusted"].rank(method="first")
        quartile[significant] = np.floor(4 * (order - 1) / n_significant) + 1
    table["significance_quartile"] = quartile
    return table


def pseudobulk_de(
    pseudobulk: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    design: str,
    contrast: tuple[str, str, str],
    lfc_null: float = 0.0,
    min_samples_per_group: int = MIN_SAMPLES_PER_GROUP,
    n_cpus: int = 1,
    alpha: float = ALPHA,
    shrink: bool = True,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Condition DE per cell type.

    Parameters
    ----------
    pseudobulk : dict[str, tuple[pd.DataFrame, pd.DataFrame]]
        Output of ``aggregate_counts``.
    design : str
        Formula, e.g. ``"~ sex + age + condition"``.
    contrast : tuple[str, str, str]
        ``(variable, tested_level, reference_level)``.
    lfc_null : float
        Minimum absolute log2 fold change tested (0 tests for any change).
    min_samples_per_group : int
        Cell types with fewer units in a compared group are not tested.
    n_cpus : int
        CPUs for pydeseq2.
    alpha : float
        FDR level (independent filtering and ``significance_quartile``).
    shrink : bool
        Whether to add shrunken log fold changes.

    Returns
    -------
    results : pl.DataFrame
        ``recipes.bulk.deseq2.RESULT_COLUMNS`` plus ``DE_EXTRA_COLUMNS``; empty with the same
        schema when nothing was testable. ``lfc_ci_low`` and ``lfc_ci_high`` are the Wald
        interval at ``CI_LEVEL`` of the unshrunken log fold change; ``rank`` orders genes by
        p-value within the cell type; ``significance_quartile`` (1-4) is the quartile of the
        adjusted p-value among the genes with ``p_adjusted < alpha`` of the cell type, missing for
        the others.
    not_testable : pl.DataFrame
        ``cell_type``, ``reason``: too few units per group, or no gene passing the design-based
        filter. Genes without a finite p-value keep a missing ``p_adjusted_global`` (BH runs over
        the finite p-values only; with none, every value is missing).
    """
    variable, tested, reference = contrast
    z = stats.norm.ppf(0.5 + CI_LEVEL / 2)
    tables, skipped = [], []
    for cell_type, (counts, meta) in pseudobulk.items():
        group_sizes = meta[variable].astype(str).value_counts()
        sizes = {level: int(group_sizes.get(level, 0)) for level in (tested, reference)}
        if min(sizes.values()) < min_samples_per_group:
            skipped.append({"cell_type": cell_type, "reason": f"units per group {sizes}"})
            continue
        genes = filter_genes(counts, meta, variable)
        if genes.size == 0:
            skipped.append({"cell_type": cell_type, "reason": "no genes pass the gene filter"})
            continue
        results, _ = run_deseq2(
            counts[genes],
            meta,
            design,
            contrast,
            lfc_null=lfc_null,
            alpha=alpha,
            n_cpus=n_cpus,
            shrink=shrink,
        )
        table = _rank_and_quartile(results.to_pandas(), alpha)
        table["lfc_ci_low"] = table["log2_fold_change"] - z * table["lfc_se"]
        table["lfc_ci_high"] = table["log2_fold_change"] + z * table["lfc_se"]
        table["cell_type"] = cell_type
        tables.append(table)
    skipped_table = pl.DataFrame(skipped, schema={"cell_type": pl.Utf8, "reason": pl.Utf8})
    columns = [*RESULT_COLUMNS, *DE_EXTRA_COLUMNS]
    if not tables:
        schema = {c: pl.Float64 for c in columns}
        schema.update({"gene_id": pl.Utf8, "contrast": pl.Utf8, "cell_type": pl.Utf8})
        return pl.DataFrame(schema=schema).select(columns), skipped_table
    results = pd.concat(tables, ignore_index=True)
    p = results["p_value"].to_numpy(dtype=float)
    p_global = np.full(p.size, np.nan)
    finite = np.isfinite(p)
    p_global[finite] = multipletests(p[finite], method="fdr_bh")[1]
    results["p_adjusted_global"] = p_global
    return pl.from_pandas(results[columns], nan_to_null=False), skipped_table


def detectability_report(
    results: pl.DataFrame,
    pseudobulk: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    genes: list[str],
    margin: float,
    contrast: tuple[str, str, str],
    alpha: float = ALPHA,
) -> pl.DataFrame:
    """Report, for genes of interest, whether each changed, did not change, or is inconclusive.

    A gene missing from a DE list is not evidence that it does not change
    (``SCRNASEQ_GUIDELINES.md`` ``SC-DE-009``). ``no_change_of_interest`` requires an equivalence
    test: two one-sided tests (TOST) that the log fold change lies within ``(-margin, margin)``,
    Holm-adjusted across the genes of the report within each cell type. ``margin`` is the smallest
    log2 fold change of interest and MUST be fixed before the results are seen.

    Parameters
    ----------
    results : pl.DataFrame
        Output of ``pseudobulk_de``.
    pseudobulk : dict[str, tuple[pd.DataFrame, pd.DataFrame]]
        The input of ``pseudobulk_de`` (to tell filtered genes and untested cell types apart).
    genes : list[str]
        Gene identifiers (``gene_id`` of the results).
    margin : float
        Equivalence margin on the log2 fold change (> 0).
    contrast : tuple[str, str, str]
        The contrast of ``pseudobulk_de``.
    alpha : float
        Level of the DE calls and of the equivalence tests.

    Returns
    -------
    pl.DataFrame
        One row per (cell type, gene): ``cell_type``, ``gene_id``, ``status`` (``changed``,
        ``no_change_of_interest``, ``inconclusive``, ``not_tested``), ``reason`` (for
        ``not_tested``: ``cell_type_not_testable``, ``gene_absent``, ``filtered_low_counts``,
        ``no_p_value``), ``base_mean``, ``log2_fold_change``, ``lfc_ci_low``, ``lfc_ci_high``,
        ``p_adjusted``, ``p_equivalence`` (Holm-adjusted TOST), ``n_units_tested``,
        ``n_units_reference`` and ``median_cells_per_unit``.

    Raises
    ------
    ValueError
        If ``margin`` is not positive.
    """
    if not margin > 0:
        raise ValueError(f"margin must be > 0, got {margin}")
    variable, tested, reference = contrast
    res = results.to_pandas().set_index(["cell_type", "gene_id"])
    rows = []
    for cell_type, (counts, meta) in pseudobulk.items():
        levels = meta[variable].astype(str)
        base = {
            "cell_type": cell_type,
            "n_units_tested": int((levels == tested).sum()),
            "n_units_reference": int((levels == reference).sum()),
            "median_cells_per_unit": float(meta["n_cells"].median()),
        }
        tested_genes = res.loc[cell_type] if cell_type in res.index.get_level_values(0) else None
        block = []
        for gene in genes:
            row = {**base, "gene_id": gene, "status": "not_tested", "reason": None}
            if tested_genes is None:
                row["reason"] = "cell_type_not_testable"
            elif gene not in counts.columns:
                row["reason"] = "gene_absent"
            elif gene not in tested_genes.index:
                row["reason"] = "filtered_low_counts"
            else:
                g = tested_genes.loc[gene]
                row.update(
                    {
                        k: float(g[k])
                        for k in (
                            "base_mean",
                            "log2_fold_change",
                            "lfc_se",
                            "lfc_ci_low",
                            "lfc_ci_high",
                            "p_value",
                            "p_adjusted",
                        )
                    }
                )
                if not np.isfinite(row["p_value"]):
                    row["reason"] = "no_p_value"
            block.append(row)
        _equivalence(block, margin, alpha)
        rows += block
    columns = {
        "cell_type": pl.Utf8,
        "gene_id": pl.Utf8,
        "status": pl.Utf8,
        "reason": pl.Utf8,
        "base_mean": pl.Float64,
        "log2_fold_change": pl.Float64,
        "lfc_ci_low": pl.Float64,
        "lfc_ci_high": pl.Float64,
        "p_adjusted": pl.Float64,
        "p_equivalence": pl.Float64,
        "n_units_tested": pl.Int64,
        "n_units_reference": pl.Int64,
        "median_cells_per_unit": pl.Float64,
    }
    return pl.DataFrame([{c: r.get(c) for c in columns} for r in rows], schema=columns)


def _equivalence(block: list[dict], margin: float, alpha: float) -> None:
    """Set ``status`` and ``p_equivalence`` of the tested genes of one cell type, in place."""
    testable = [r for r in block if r["reason"] is None]
    if not testable:
        return
    lfc = np.array([r["log2_fold_change"] for r in testable])
    se = np.array([r["lfc_se"] for r in testable])
    p_lower = stats.norm.sf((lfc + margin) / se)  # H0: lfc <= -margin
    p_upper = stats.norm.sf((margin - lfc) / se)  # H0: lfc >= margin
    p_tost = multipletests(np.maximum(p_lower, p_upper), method="holm")[1]
    for row, p_eq in zip(testable, p_tost, strict=True):
        row["p_equivalence"] = float(p_eq)
        p_adj = row["p_adjusted"]
        if np.isfinite(p_adj) and p_adj < alpha:
            row["status"] = "changed"
        elif p_eq < alpha:
            row["status"] = "no_change_of_interest"
        else:
            row["status"] = "inconclusive"


def null_split_de(
    pseudobulk: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    design: str,
    contrast: tuple[str, str, str],
    mode: str = "split_controls",
    n_splits: int = N_NULL_SPLITS,
    strata_key: str | None = None,
    null_design: str | None = None,
    min_samples_per_group: int = MIN_SAMPLES_PER_GROUP,
    alpha: float = ALPHA,
    seed: int = 0,
    n_cpus: int = 1,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Global-null diagnostic of the pseudobulk DE pipeline, per cell type (``SC-DE-007``).

    ``mode="split_controls"`` keeps the units of the reference level of the contrast and splits
    them at random into two groups of equal size (``NULL_SPLIT_COLUMN``), re-running the whole DE
    each time. ``mode="permute_labels"`` permutes the condition labels of all units, as
    ``RNASEQ_BULK_GUIDELINES.md`` ``BULK-DE-005``. In both, every discovery is false; passing shows
    that the pipeline does not create discoveries from nothing, not that it controls the FDR. The
    caller is responsible for exchangeability (``STATISTICS_GUIDELINES.md`` §10).

    Parameters
    ----------
    pseudobulk : dict[str, tuple[pd.DataFrame, pd.DataFrame]]
        Output of ``aggregate_counts``.
    design : str
        Design of the real analysis.
    contrast : tuple[str, str, str]
        Contrast of the real analysis; its reference level defines the controls.
    mode : str
        ``"split_controls"`` or ``"permute_labels"``.
    n_splits : int
        Random splits (or permutations) per cell type.
    strata_key : str | None
        Unit-level column within which labels are permuted (batch); None permutes freely.
    null_design : str | None
        Design for ``split_controls``; by default ``design`` with the contrast variable replaced
        by ``NULL_SPLIT_COLUMN``. Give it when a covariate is constant among the controls.
    min_samples_per_group : int
        Cell types whose groups would have fewer units are reported as not testable.
    alpha : float
        FDR level of the real analysis.
    seed : int
        Seed of the splits.
    n_cpus : int
        CPUs for pydeseq2.

    Returns
    -------
    per_split : pl.DataFrame
        ``cell_type`` plus the per-permutation table of ``recipes.bulk.fdr_calibration``.
    summary : pl.DataFrame
        One row per cell type: ``cell_type``, ``mode``, ``n_units``, ``status`` (``tested`` or the
        reason it was not), and the summary of ``fdr_calibration``.

    Raises
    ------
    ValueError
        If ``mode`` is unknown.
    """
    if mode not in ("split_controls", "permute_labels"):
        raise ValueError(f"mode must be 'split_controls' or 'permute_labels', got {mode!r}")
    variable, tested, reference = contrast
    tables, summaries = [], []
    for index, (cell_type, (counts, meta)) in enumerate(pseudobulk.items()):
        if mode == "split_controls":
            units = meta.index[meta[variable].astype(str) == reference]
            null_meta = meta.loc[units].copy()
            null_meta[NULL_SPLIT_COLUMN] = [
                NULL_SPLIT_LEVELS[i % 2] for i in range(len(null_meta))
            ]
            column, level_b, level_a = NULL_SPLIT_COLUMN, *NULL_SPLIT_LEVELS[::-1]
            formula = null_design or re.sub(
                rf"\b{re.escape(variable)}\b", NULL_SPLIT_COLUMN, design
            )
            group_size = len(null_meta) // 2
        else:
            null_meta = meta.copy()
            column, level_b, level_a, formula = variable, tested, reference, design
            group_size = int(null_meta[variable].astype(str).value_counts().min())
        summary = {"cell_type": cell_type, "mode": mode, "n_units": len(null_meta)}
        if group_size < min_samples_per_group:
            summaries.append({**summary, "status": f"not_testable: {group_size} units per group"})
            continue
        null_counts = counts.loc[null_meta.index]

        def pipeline(m, c=null_counts, col=column, f=formula, lv=(level_b, level_a)):
            genes = filter_genes(c, m, col)
            res, _ = run_deseq2(
                c[genes], m, f, (col, *lv), alpha=alpha, n_cpus=n_cpus, shrink=False
            )
            return res["p_adjusted"].to_numpy()

        table, result = fdr_calibration(
            pipeline,
            null_meta,
            column,
            strata_col=strata_key,
            n_permutations=n_splits,
            alpha=alpha,
            seed=seed + index,
        )
        tables.append(table.with_columns(cell_type=pl.lit(cell_type)))
        summaries.append({**summary, "status": "tested", **result})
    per_split = (
        pl.concat(tables).select("cell_type", pl.exclude("cell_type"))
        if tables
        else pl.DataFrame(schema={"cell_type": pl.Utf8})
    )
    return per_split, pl.DataFrame(summaries)


def min_cells_diagnostic(
    pseudobulk: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    n_pcs: int = DIAGNOSTIC_PCS,
    n_genes: int = DIAGNOSTIC_GENES,
    bins: tuple[int, ...] = CELL_BINS,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Principal components of pseudobulk samples against their number of cells (``SC-DE-002``).

    Run on an aggregation with ``min_cells=1``, before any DE result is seen. Per cell type, log
    CPM of the ``n_genes`` most variable genes (mean CPM at least ``DIAGNOSTIC_MIN_CPM``), centred;
    the first ``n_pcs`` components; and the distance of each unit from the median unit, in robust
    z-scores (median and scaled MAD of each component). Samples with few cells that drift away from
    the rest show as a larger median distance in the low bins and a negative rank correlation
    between distance and number of cells. The threshold is a judgement fixed from this table; the
    function does not choose it.

    Parameters
    ----------
    pseudobulk : dict[str, tuple[pd.DataFrame, pd.DataFrame]]
        Output of ``aggregate_counts`` (``min_cells=1``).
    n_pcs : int
        Components used.
    n_genes : int
        Most variable genes used.
    bins : tuple[int, ...]
        Lower edges of the cells-per-unit bins.

    Returns
    -------
    units : pl.DataFrame
        ``cell_type``, ``unit``, ``n_cells``, ``total_counts``, ``pc1``, ..., ``distance``.
    by_bin : pl.DataFrame
        ``cell_type``, ``cells_bin`` (e.g. ``"10-19"``), ``n_units``, ``median_distance``, and per
        cell type ``rho_distance_log_cells`` and ``rho_pc1_log_cells`` (Spearman).
    """
    unit_rows, bin_rows = [], []
    edges = list(bins)
    labels = [
        f"{lo}-{hi - 1}" if hi - 1 > lo else f"{lo}"
        for lo, hi in zip(edges[:-1], edges[1:], strict=True)
    ] + [f"{edges[-1]}+"]
    for cell_type, (counts, meta) in pseudobulk.items():
        values = counts.to_numpy(dtype=float)
        total = values.sum(axis=1)
        if values.shape[0] < 3 or np.any(total == 0):
            logger.warning(
                f"{cell_type}: not evaluated ({values.shape[0]} units, "
                f"{int(np.sum(total == 0))} with no counts)"
            )
            continue
        cpm = values / total[:, None] * 1e6
        keep = cpm.mean(axis=0) >= DIAGNOSTIC_MIN_CPM
        log_cpm = np.log2(cpm[:, keep] + 1)
        top = np.argsort(log_cpm.var(axis=0))[::-1][:n_genes]
        centred = log_cpm[:, top] - log_cpm[:, top].mean(axis=0)
        u, s, _ = np.linalg.svd(centred, full_matrices=False)
        k = min(n_pcs, s.size)
        pcs = u[:, :k] * s[:k]
        med = np.median(pcs, axis=0)
        mad = stats.median_abs_deviation(pcs, axis=0, scale="normal")
        mad[mad == 0] = 1.0
        distance = np.sqrt((((pcs - med) / mad) ** 2).sum(axis=1))
        n_cells = meta["n_cells"].to_numpy()
        for i, unit in enumerate(counts.index):
            unit_rows.append(
                {
                    "cell_type": cell_type,
                    "unit": str(unit),
                    "n_cells": int(n_cells[i]),
                    "total_counts": float(total[i]),
                    **{f"pc{j + 1}": float(pcs[i, j]) for j in range(k)},
                    "distance": float(distance[i]),
                }
            )
        log_cells = np.log(n_cells)
        rho_distance = stats.spearmanr(log_cells, distance).statistic
        rho_pc1 = stats.spearmanr(log_cells, pcs[:, 0]).statistic
        which = np.searchsorted(edges, n_cells, side="right") - 1
        for b, label in enumerate(labels):
            in_bin = which == b
            bin_rows.append(
                {
                    "cell_type": cell_type,
                    "cells_bin": label,
                    "n_units": int(in_bin.sum()),
                    "median_distance": float(np.median(distance[in_bin]))
                    if in_bin.any()
                    else None,
                    "rho_distance_log_cells": float(rho_distance),
                    "rho_pc1_log_cells": float(rho_pc1),
                }
            )
    return pl.DataFrame(unit_rows), pl.DataFrame(bin_rows)


def power_inputs(
    pseudobulk: dict[str, tuple[pd.DataFrame, pd.DataFrame]], variable: str
) -> pl.DataFrame:
    """Units, cells per unit and UMIs per cell, per cell type and group (``SC-DE-008``).

    Parameters
    ----------
    pseudobulk : dict[str, tuple[pd.DataFrame, pd.DataFrame]]
        Output of ``aggregate_counts``.
    variable : str
        Unit-level column with the compared groups.

    Returns
    -------
    pl.DataFrame
        ``cell_type``, ``group``, ``n_units``, ``cells_per_unit_median``, ``cells_per_unit_q1``,
        ``cells_per_unit_q3``, ``umis_per_cell_median`` (median over units of the mean UMIs per
        cell) and ``umis_per_unit_median``.
    """
    rows = []
    for cell_type, (counts, meta) in pseudobulk.items():
        total = counts.to_numpy(dtype=float).sum(axis=1)
        n_cells = meta["n_cells"].to_numpy(dtype=float)
        groups = meta[variable].astype(str).to_numpy()
        for group in sorted(set(groups)):
            g = groups == group
            rows.append(
                {
                    "cell_type": cell_type,
                    "group": group,
                    "n_units": int(g.sum()),
                    "cells_per_unit_median": float(np.median(n_cells[g])),
                    "cells_per_unit_q1": float(np.quantile(n_cells[g], 0.25)),
                    "cells_per_unit_q3": float(np.quantile(n_cells[g], 0.75)),
                    "umis_per_cell_median": float(np.median(total[g] / n_cells[g])),
                    "umis_per_unit_median": float(np.median(total[g])),
                }
            )
    return pl.DataFrame(rows)
