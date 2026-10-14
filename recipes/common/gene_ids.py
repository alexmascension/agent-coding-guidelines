"""Recipe: gene annotation table from the GTF, stable Ensembl IDs, and marker dictionaries resolved to IDs.

Guidelines: ``SCRNASEQ_GUIDELINES.md`` §1.3 (gene identifiers and symbols) and §4.2 (mitochondrial genes by
chromosome); ``RNASEQ_BULK_GUIDELINES.md`` §2 (reference genome and annotation).

Procedure
---------
1. ``annotation_table_from_gtf``: one row per gene from the GTF used for quantification (``gene_id`` stable, without
   version; ``gene_id_version``; ``gene_symbol``; ``chromosome``; ``gene_biotype``). Works for Ensembl GTFs (version in
   the ``gene_version`` attribute) and GENCODE GTFs (version inside ``gene_id``).
2. ``set_var_gene_ids``: re-index ``adata.var`` by stable ID and add ``gene_id_version`` and ``gene_symbol``.
3. ``resolve_markers``: translate a marker dictionary written in symbols into IDs, with a report of every symbol that
   was not found, matched only with a different capitalization, or matched several IDs.

Adapting
--------
* The rule for symbols with several IDs (``multi_id_rule``) is a project decision recorded in ``PROJECT_CONTEXT.md``.
* Case-insensitive matches are reported, never used: human symbols are not mouse symbols (§1.3).
"""

import gzip
import logging
import re
from pathlib import Path
from typing import Literal

import anndata as ad
import numpy as np
import polars as pl

logger = logging.getLogger(__name__)

MITOCHONDRIAL_CHROMOSOMES = ("MT", "chrM", "M")  # Ensembl, UCSC/GENCODE, some custom references
_ATTRIBUTE = re.compile(r'(\S+) "([^"]*)"')
_VERSION = re.compile(r"\.\d+")  # first ".N"; keeps GENCODE suffixes such as "_PAR_Y"

ANNOTATION_SCHEMA = {
    "gene_id": pl.Utf8,
    "gene_id_version": pl.Utf8,
    "gene_symbol": pl.Utf8,
    "chromosome": pl.Utf8,
    "gene_biotype": pl.Utf8,
}
REPORT_COLUMNS = ("cell_type", "symbol", "status", "gene_ids", "candidates")


def strip_version(gene_id: str) -> str:
    """Remove the version from an Ensembl ID (``ENSG00000141510.18`` -> ``ENSG00000141510``)."""
    return _VERSION.sub("", gene_id, count=1)


def annotation_table_from_gtf(gtf_path: str | Path) -> pl.DataFrame:
    """Read gene records of a GTF (plain or gzip) into an annotation table.

    Parameters
    ----------
    gtf_path : str | Path
        GTF used for quantification (the same file, so that IDs and versions match the counts).

    Returns
    -------
    pl.DataFrame
        Columns of ``ANNOTATION_SCHEMA``, one row per gene. Missing symbols are null, never replaced by the ID.

    Raises
    ------
    ValueError
        If stable IDs are not unique after removing versions.
    """
    path = Path(gtf_path)
    opener = gzip.open if path.suffix == ".gz" else open
    rows = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = dict(_ATTRIBUTE.findall(fields[8]))
            raw_id = attributes["gene_id"]
            if "gene_version" in attributes and "." not in raw_id:
                versioned = f"{raw_id}.{attributes['gene_version']}"
            else:
                versioned = raw_id
            rows.append(
                {
                    "gene_id": strip_version(raw_id),
                    "gene_id_version": versioned,
                    "gene_symbol": attributes.get("gene_name"),
                    "chromosome": fields[0],
                    "gene_biotype": attributes.get("gene_biotype", attributes.get("gene_type")),
                }
            )
    table = pl.DataFrame(rows, schema=ANNOTATION_SCHEMA)
    duplicated = table.filter(pl.col("gene_id").is_duplicated())["gene_id"].unique().to_list()
    if duplicated:
        raise ValueError(f"stable IDs not unique after removing versions: {duplicated[:5]}")
    logger.info(f"annotation table: {table.height} genes from {path.name}")
    return table


def set_var_gene_ids(adata: ad.AnnData, annotation: pl.DataFrame) -> pl.DataFrame:
    """Index ``adata.var`` by stable Ensembl ID and add ``gene_id_version`` and ``gene_symbol`` (in place).

    The current index may contain versioned or unversioned IDs. Genes absent from the annotation keep null symbol and
    version and are returned in the report.

    Parameters
    ----------
    adata : ad.AnnData
        Object whose ``var_names`` are Ensembl IDs.
    annotation : pl.DataFrame
        Output of ``annotation_table_from_gtf``.

    Returns
    -------
    pl.DataFrame
        IDs of ``adata`` not found in the annotation (column ``gene_id``); empty if all matched.

    Raises
    ------
    ValueError
        If ``var_names`` are not unique after removing versions.
    """
    stable = np.array([strip_version(g) for g in adata.var_names])
    if len(set(stable)) != stable.size:
        raise ValueError("var_names are not unique after removing versions")
    lookup = annotation.select("gene_id", "gene_id_version", "gene_symbol").to_pandas().set_index("gene_id")
    matched = lookup.reindex(stable)
    adata.var_names = stable
    adata.var["gene_id_version"] = matched["gene_id_version"].to_numpy()
    adata.var["gene_symbol"] = matched["gene_symbol"].to_numpy()
    absent = stable[matched["gene_id_version"].isna().to_numpy()]
    missing = pl.DataFrame({"gene_id": absent}, schema={"gene_id": pl.Utf8})
    if missing.height:
        logger.warning(f"{missing.height} genes of the matrix are absent from the annotation")
    return missing


def resolve_markers(
    markers: dict[str, list[str]],
    annotation: pl.DataFrame,
    multi_id_rule: Literal["keep_all", "most_expressed"] = "keep_all",
    expression: dict[str, float] | None = None,
) -> tuple[dict[str, list[str]], pl.DataFrame]:
    """Translate a marker dictionary from symbols to stable Ensembl IDs and report every problem.

    Parameters
    ----------
    markers : dict[str, list[str]]
        Cell type -> marker symbols (species-specific).
    annotation : pl.DataFrame
        Output of ``annotation_table_from_gtf`` (same release as the counts).
    multi_id_rule : {"keep_all", "most_expressed"}
        What to do with a symbol that maps to several IDs (project decision).
    expression : dict[str, float] | None
        Gene ID -> mean expression; required for ``"most_expressed"``.

    Returns
    -------
    resolved : dict[str, list[str]]
        Cell type -> marker IDs (order preserved, no duplicates).
    report : pl.DataFrame
        One row per (cell type, symbol): ``cell_type``, ``symbol``, ``status`` (``resolved``, ``multiple_ids``,
        ``case_mismatch``, ``not_found``), ``gene_ids`` (semicolon-separated IDs kept), ``candidates`` (all IDs, or
        the symbols that match only with different capitalization).

    Raises
    ------
    ValueError
        If ``multi_id_rule="most_expressed"`` and ``expression`` is missing.
    """
    if multi_id_rule == "most_expressed" and expression is None:
        raise ValueError("multi_id_rule='most_expressed' needs the expression mapping")
    by_symbol: dict[str, list[str]] = {}
    for gene_id, symbol in annotation.select("gene_id", "gene_symbol").drop_nulls().iter_rows():
        by_symbol.setdefault(symbol, []).append(gene_id)
    by_upper: dict[str, list[str]] = {}
    for symbol in by_symbol:
        by_upper.setdefault(symbol.upper(), []).append(symbol)

    resolved: dict[str, list[str]] = {}
    report_rows = []
    for cell_type, symbols in markers.items():
        kept: list[str] = []
        for symbol in symbols:
            ids = by_symbol.get(symbol, [])
            if len(ids) == 1:
                status, chosen, candidates = "resolved", ids, ids
            elif len(ids) > 1:
                status, candidates = "multiple_ids", ids
                if multi_id_rule == "keep_all":
                    chosen = ids
                else:
                    chosen = [max(ids, key=lambda g: expression.get(g, float("-inf")))]
            else:
                chosen = []
                near = [s for s in by_upper.get(symbol.upper(), []) if s != symbol]
                status, candidates = ("case_mismatch", near) if near else ("not_found", [])
            kept.extend(g for g in chosen if g not in kept)
            report_rows.append(
                {
                    "cell_type": cell_type,
                    "symbol": symbol,
                    "status": status,
                    "gene_ids": ";".join(chosen),
                    "candidates": ";".join(candidates),
                }
            )
        resolved[cell_type] = kept
    report = pl.DataFrame(report_rows, schema=dict.fromkeys(REPORT_COLUMNS, pl.Utf8))
    problems = report.filter(pl.col("status") != "resolved")
    if problems.height:
        logger.warning(f"{problems.height} marker symbols not resolved one-to-one; see the report")
    return resolved, report


def symbols_for(gene_ids: list[str], annotation: pl.DataFrame) -> list[str]:
    """Return display labels for IDs: the symbol when it exists, otherwise the ID (for plots and tables, §1.3)."""
    lookup = dict(annotation.select("gene_id", "gene_symbol").iter_rows())
    return [lookup.get(g) or g for g in gene_ids]
