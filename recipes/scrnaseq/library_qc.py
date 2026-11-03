"""Recipe: library-level QC from STARsolo summaries, compared across the libraries of a project.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §3 (``SC-QUANT-006``, ``SC-QUANT-007``) and §1.2
(``SC-DATA-011``, ``uns["library_qc"]``).

Procedure
---------
1. ``nfcore_scrnaseq_paths``: locate ``Summary.csv``, ``Barcodes.stats`` and ``Log.out`` of one
   sample in the output layout of ``nf-core/scrnaseq``
   (``<outdir>/star/<sample>/<sample>.Solo.out/...``).
2. ``read_summary_csv``: the metrics of ``Solo.out/<feature>/Summary.csv``. STARsolo writes the
   feature name into some row names (``Median GeneFull_Ex50pAS per Cell``); the feature is
   detected and the rows are mapped to fixed column names. ``NoMulti`` (multi-gene reads not
   counted) becomes a missing value.
3. ``read_barcodes_stats``: the counts of ``Solo.out/Barcodes.stats`` (STAR 2.7.11b names, e.g.
   ``noNoWLmatch``, ``yesWLmatchExact``); the fraction of reads whose barcode matches no
   whitelist entry points to a wrong chemistry or whitelist.
4. ``library_qc_table``: one row per ``sample_id`` with the metrics, the number of cells loaded
   when given, the source files, the feature and the STAR version.
5. ``flag_libraries``: per metric, flags libraries more than ``n_mads`` scaled MADs from the
   median of the libraries, in the direction that indicates a problem (e.g. a low fraction of
   valid barcodes).
6. ``store_library_qc``: writes the table to ``adata.uns["library_qc"]`` after the captures are
   concatenated, checking that every ``sample_id`` of ``obs`` has a row.

Policies that the caller sees in the report
-------------------------------------------
* Fewer than ``MIN_LIBRARIES`` libraries with a value: the metric is ``not_evaluated_too_few``.
* MAD = 0 across libraries: ``not_evaluated_mad_zero``; no library is flagged on that metric.
* A metric absent from every summary: ``not_available``; the column stays missing, never 0.
* Flags are ``None`` where the metric was not evaluated, so a missing comparison is not read as
  a pass.

Defaults
--------
``N_MADS_LIBRARY = 3`` scaled MADs and ``MIN_LIBRARIES = 3`` are project defaults: they compare
libraries with each other and do not replace the judgement of ``SC-QUANT-006`` (investigate,
then decide).
"""

import logging
import re
from collections.abc import Mapping
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl

from recipes.scrnaseq.qc import MAD_SCALE

logger = logging.getLogger(__name__)

N_MADS_LIBRARY = 3
MIN_LIBRARIES = 3
UNS_KEY = "library_qc"

# Summary.csv row -> column. "{f}" stands for the feature name STARsolo writes into the row.
SUMMARY_ROWS = {
    "Number of Reads": "n_reads",
    "Reads With Valid Barcodes": "frac_valid_barcodes",
    "Sequencing Saturation": "sequencing_saturation",
    "Q30 Bases in CB+UMI": "q30_cb_umi",
    "Q30 Bases in RNA read": "q30_rna",
    "Reads Mapped to Genome: Unique+Multiple": "frac_mapped_genome",
    "Reads Mapped to Genome: Unique": "frac_mapped_genome_unique",
    "Reads Mapped to {f}: Unique+Multiple {f}": "frac_mapped_feature",
    "Reads Mapped to {f}: Unique {f}": "frac_mapped_feature_unique",
    "Estimated Number of Cells": "estimated_cells",
    "Unique Reads in Cells Mapped to {f}": "unique_reads_in_cells",
    "Fraction of Unique Reads in Cells": "frac_unique_reads_in_cells",
    "Mean Reads per Cell": "mean_reads_per_cell",
    "Median Reads per Cell": "median_reads_per_cell",
    "UMIs in Cells": "umis_in_cells",
    "Mean UMI per Cell": "mean_umi_per_cell",
    "Median UMI per Cell": "median_umi_per_cell",
    "Mean {f} per Cell": "mean_genes_per_cell",
    "Median {f} per Cell": "median_genes_per_cell",
    "Total {f} Detected": "total_genes_detected",
}
SUMMARY_COLUMNS = list(SUMMARY_ROWS.values())

# Metrics compared across libraries and the direction that indicates a problem.
FLAG_DIRECTIONS = {
    "n_reads": "both",
    "frac_valid_barcodes": "lower",
    "frac_cb_no_whitelist_match": "upper",
    "sequencing_saturation": "lower",
    "q30_cb_umi": "lower",
    "q30_rna": "lower",
    "frac_mapped_genome": "lower",
    "frac_mapped_feature_unique": "lower",
    "estimated_cells": "both",
    "cells_over_loaded": "both",
    "frac_unique_reads_in_cells": "lower",
    "median_umi_per_cell": "both",
    "median_genes_per_cell": "both",
}

REPORT_SCHEMA = {
    "metric": pl.Utf8,
    "direction": pl.Utf8,
    "n_libraries": pl.Int64,
    "median": pl.Float64,
    "scaled_mad": pl.Float64,
    "lower": pl.Float64,
    "upper": pl.Float64,
    "status": pl.Utf8,
    "n_flagged": pl.Int64,
}

_FEATURE_ROW = re.compile(r"^Reads Mapped to (?!Genome:)(\S+): Unique (\S+)$")
_STAR_VERSION = re.compile(r"STAR version[=\s]+(\S+)")


def nfcore_scrnaseq_paths(
    outdir: str | Path, sample: str, feature: str = "Gene"
) -> dict[str, Path]:
    """Return the STARsolo files of one sample in the ``nf-core/scrnaseq`` output layout.

    Parameters
    ----------
    outdir : str | Path
        The ``--outdir`` of the pipeline run (``--aligner star``).
    sample : str
        Sample identifier as in the samplesheet (the pipeline's ``meta.id``).
    feature : str
        Gene-level feature whose ``Summary.csv`` is read (``Gene`` or ``GeneFull_Ex50pAS``).

    Returns
    -------
    dict[str, Path]
        ``summary``, ``barcodes_stats`` and ``log_out`` paths (not checked for existence).
    """
    sample_dir = Path(outdir) / "star" / sample
    solo = sample_dir / f"{sample}.Solo.out"
    return {
        "summary": solo / feature / "Summary.csv",
        "barcodes_stats": solo / "Barcodes.stats",
        "log_out": sample_dir / f"{sample}.Log.out",
    }


def _to_number(text: str) -> float | None:
    try:
        return float(text)
    except ValueError:
        return None


def read_summary_csv(path: str | Path) -> dict[str, float | str | None]:
    """Read a STARsolo ``Summary.csv`` into fixed column names.

    Parameters
    ----------
    path : str | Path
        ``Solo.out/<feature>/Summary.csv``.

    Returns
    -------
    dict[str, float | str | None]
        One entry per column of ``SUMMARY_COLUMNS`` (``None`` when the row is absent or not
        numeric, e.g. ``NoMulti``) and ``feature``, the feature name found in the file.

    Raises
    ------
    ValueError
        If a line does not have the ``name,value`` form.
    """
    rows: dict[str, str] = {}
    for number, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        name, sep, value = line.rpartition(",")
        if not sep:
            raise ValueError(f"{path}:{number}: expected 'name,value', got {line!r}")
        rows[name.strip()] = value.strip()

    feature = None
    for name in rows:
        match = _FEATURE_ROW.match(name)
        if match and match.group(1) == match.group(2):
            feature = match.group(1)
            break
    if feature is None:
        logger.warning(f"{path}: feature name not found; feature-specific metrics are missing")

    values: dict[str, float | str | None] = {"feature": feature}
    for template, column in SUMMARY_ROWS.items():
        if "{f}" in template and feature is None:
            values[column] = None
            continue
        raw = rows.get(template.format(f=feature))
        values[column] = None if raw is None else _to_number(raw)
    missing = [column for column in SUMMARY_COLUMNS if values[column] is None]
    if missing:
        logger.info(f"{path}: no numeric value for {missing}")
    return values


def read_barcodes_stats(path: str | Path) -> dict[str, int]:
    """Read ``Solo.out/Barcodes.stats`` (one ``name count`` pair per line).

    Parameters
    ----------
    path : str | Path
        ``Solo.out/Barcodes.stats``.

    Returns
    -------
    dict[str, int]
        Count per category as written by STAR.

    Raises
    ------
    ValueError
        If a line does not have two fields with an integer count.
    """
    counts: dict[str, int] = {}
    for number, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), start=1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 2 or not fields[1].isdigit():
            raise ValueError(f"{path}:{number}: expected 'name count', got {line!r}")
        counts[fields[0]] = int(fields[1])
    return counts


def star_version(log_out: str | Path) -> str | None:
    """Return the STAR version recorded in ``Log.out``, or ``None`` if absent."""
    path = Path(log_out)
    if not path.is_file():
        return None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for _, line in zip(range(50), handle, strict=False):
            match = _STAR_VERSION.search(line)
            if match:
                return match.group(1)
    return None


def library_qc_table(
    samples: Mapping[str, Mapping[str, str | Path]],
    cells_loaded: Mapping[str, float] | None = None,
) -> pl.DataFrame:
    """Build the library QC table, one row per sample.

    Parameters
    ----------
    samples : Mapping[str, Mapping[str, str | Path]]
        ``sample_id`` -> paths with keys ``summary`` and ``barcodes_stats`` (required) and
        ``log_out`` (optional), e.g. from ``nfcore_scrnaseq_paths``.
    cells_loaded : Mapping[str, float] | None
        Cells loaded per ``sample_id``, when known; gives ``cells_over_loaded``.

    Returns
    -------
    pl.DataFrame
        Metrics of ``SUMMARY_COLUMNS``, ``frac_cb_no_whitelist_match``, ``cells_loaded``,
        ``cells_over_loaded``, ``feature``, ``star_version`` and the source files, sorted by
        ``sample_id``.

    Raises
    ------
    FileNotFoundError
        If a required file is missing.
    """
    records = []
    for sample_id in sorted(samples):
        paths = samples[sample_id]
        for key in ("summary", "barcodes_stats"):
            if not Path(paths[key]).is_file():
                raise FileNotFoundError(f"{sample_id}: {key} file not found: {paths[key]}")
        record: dict[str, object] = {"sample_id": sample_id}
        record.update(read_summary_csv(paths["summary"]))
        barcodes = read_barcodes_stats(paths["barcodes_stats"])
        total = sum(barcodes.values())
        no_match = barcodes.get("noNoWLmatch")
        record["frac_cb_no_whitelist_match"] = (
            no_match / total if no_match is not None and total > 0 else None
        )
        loaded = None if cells_loaded is None else cells_loaded.get(sample_id)
        record["cells_loaded"] = None if loaded is None else float(loaded)
        cells = record["estimated_cells"]
        record["cells_over_loaded"] = (
            cells / loaded if cells is not None and loaded is not None and loaded > 0 else None
        )
        log_out = paths.get("log_out")
        record["star_version"] = None if log_out is None else star_version(log_out)
        record["summary_file"] = str(paths["summary"])
        record["barcodes_file"] = str(paths["barcodes_stats"])
        records.append(record)

    schema = {"sample_id": pl.Utf8, "feature": pl.Utf8}
    schema |= {column: pl.Float64 for column in SUMMARY_COLUMNS}
    schema |= {
        "frac_cb_no_whitelist_match": pl.Float64,
        "cells_loaded": pl.Float64,
        "cells_over_loaded": pl.Float64,
        "star_version": pl.Utf8,
        "summary_file": pl.Utf8,
        "barcodes_file": pl.Utf8,
    }
    return pl.DataFrame(records, schema=schema)


def flag_libraries(
    table: pl.DataFrame,
    n_mads: float = N_MADS_LIBRARY,
    min_libraries: int = MIN_LIBRARIES,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Flag libraries that depart from the others on each metric of ``FLAG_DIRECTIONS``.

    Parameters
    ----------
    table : pl.DataFrame
        Output of ``library_qc_table``.
    n_mads : float
        Width of the accepted interval in scaled MADs (MAD x 1.4826) around the median of the
        libraries.
    min_libraries : int
        Minimum number of libraries with a value for a metric to be evaluated.

    Returns
    -------
    tuple[pl.DataFrame, pl.DataFrame]
        The table with one ``flag_<metric>`` column per metric (``None`` where not evaluated) and
        ``library_flagged`` (any flag true), and a report with one row per metric.
    """
    flagged = table.clone()
    report = []
    for metric, direction in FLAG_DIRECTIONS.items():
        values = flagged.get_column(metric).to_numpy().astype(float)
        finite = np.isfinite(values)
        n_libraries = int(finite.sum())
        median = mad = lower = upper = None
        flags: list[bool | None] = [None] * flagged.height
        if n_libraries == 0:
            status = "not_available"
        elif n_libraries < min_libraries:
            status = "not_evaluated_too_few"
        else:
            median = float(np.median(values[finite]))
            mad = MAD_SCALE * float(np.median(np.abs(values[finite] - median)))
            if not mad > 0:
                status = "not_evaluated_mad_zero"
            else:
                status = "evaluated"
                lower = median - n_mads * mad if direction in ("lower", "both") else None
                upper = median + n_mads * mad if direction in ("upper", "both") else None
                for index, value in enumerate(values):
                    if not np.isfinite(value):
                        continue
                    low = lower is not None and value < lower
                    high = upper is not None and value > upper
                    flags[index] = bool(low or high)
        flagged = flagged.with_columns(pl.Series(f"flag_{metric}", flags, dtype=pl.Boolean))
        n_flagged = sum(1 for flag in flags if flag)
        if n_flagged:
            logger.warning(f"{n_flagged} libraries depart from the others on {metric}")
        report.append(
            {
                "metric": metric,
                "direction": direction,
                "n_libraries": n_libraries,
                "median": median,
                "scaled_mad": mad,
                "lower": lower,
                "upper": upper,
                "status": status,
                "n_flagged": n_flagged,
            }
        )
    flag_columns = [f"flag_{metric}" for metric in FLAG_DIRECTIONS]
    flagged = flagged.with_columns(
        pl.any_horizontal(pl.col(flag_columns).fill_null(False)).alias("library_flagged")
    )
    return flagged, pl.DataFrame(report, schema=REPORT_SCHEMA)


def store_library_qc(
    adata: ad.AnnData, table: pl.DataFrame, sample_key: str = "sample_id"
) -> None:
    """Store the library table in ``adata.uns["library_qc"]`` (in place).

    Call it after the captures are concatenated: ``anndata.concat`` keeps no ``uns`` entry by
    default.

    Parameters
    ----------
    adata : ad.AnnData
        Concatenated object with ``obs[sample_key]``.
    table : pl.DataFrame
        Output of ``library_qc_table`` or ``flag_libraries``.
    sample_key : str
        Column of ``obs`` that holds the sample identifier.

    Raises
    ------
    ValueError
        If a sample of ``obs`` has no row in the table, or the table has duplicated samples.
    """
    if sample_key not in adata.obs:
        raise ValueError(f"obs has no column {sample_key!r}")
    samples_table = table.get_column("sample_id").to_list()
    if len(set(samples_table)) != len(samples_table):
        raise ValueError("library table has duplicated sample_id values")
    samples_obs = set(adata.obs[sample_key].astype(str).unique())
    missing = sorted(samples_obs - set(samples_table))
    if missing:
        raise ValueError(f"samples in obs without library metrics: {missing}")
    extra = sorted(set(samples_table) - samples_obs)
    if extra:
        logger.warning(f"library metrics for samples absent from obs: {extra}")
    frame = _to_pandas(table.drop("sample_id"))
    frame.index = pd.Index(
        [str(sample) for sample in samples_table], dtype=object, name="sample_id"
    )
    adata.uns[UNS_KEY] = frame


def _to_pandas(table: pl.DataFrame) -> pd.DataFrame:
    """Convert to pandas with dtypes that ``write_h5ad`` stores with missing values.

    Booleans become nullable booleans; strings become categoricals, because nullable string arrays
    are refused by default in anndata versions still in use.
    """
    frame = table.to_pandas()
    for name, dtype in table.schema.items():
        values = table.get_column(name).to_list()
        if dtype == pl.Boolean:
            frame[name] = pd.array(values, dtype="boolean")
        elif dtype == pl.Utf8:
            frame[name] = pd.Categorical(values)
    return frame
