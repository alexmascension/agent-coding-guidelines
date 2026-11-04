"""Recipe: default ribosomal and haemoglobin gene lists shipped with the recipes.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §4.2 (``SC-QC-006``).

The lists are tab-separated files in ``recipes/common/data/gene_groups/`` named
``<species>_<group>.tsv`` (``human`` or ``mouse``; ``ribosomal`` or ``haemoglobin``), with columns
``gene_id`` (stable Ensembl ID without version) and ``symbol``. Header lines start with ``#`` and
hold ``key: value`` metadata: source, download date, checksum of the downloaded file and annotation
release. ``build_gene_groups.py`` writes them; they are not typed by hand.

Contents (``SC-QC-006``): cytoplasmic ribosomal proteins of the RPL and RPS families, without
mitochondrial ribosomal proteins (MRPL, MRPS) or the RPS6K kinases; haemoglobin subunits, without
HBEGF or HBS1L. ``excluded_symbols`` finds symbols that break that definition.
"""

import re
from importlib import resources
from pathlib import Path

import polars as pl

SPECIES = ("human", "mouse")
GROUPS = ("ribosomal", "haemoglobin")
# Symbols that prefix rules capture but that are not in the groups (SC-QC-006).
EXCLUDED_SYMBOL_PATTERNS = (r"^RPS6K", r"^MRP[LS]\d", r"^HBEGF$", r"^HBS1L$")
_EXCLUDED = [re.compile(pattern, re.IGNORECASE) for pattern in EXCLUDED_SYMBOL_PATTERNS]


def gene_group_path(species: str, group: str) -> Path:
    """Return the path of a shipped list (not checked for existence).

    Raises
    ------
    ValueError
        If ``species`` or ``group`` is not one of ``SPECIES`` or ``GROUPS``.
    """
    if species not in SPECIES:
        raise ValueError(f"species must be one of {SPECIES}, got {species!r}")
    if group not in GROUPS:
        raise ValueError(f"group must be one of {GROUPS}, got {group!r}")
    base = resources.files("recipes.common").joinpath("data", "gene_groups")
    return Path(str(base)) / f"{species}_{group}.tsv"


def read_gene_group(path: str | Path) -> tuple[pl.DataFrame, dict[str, str]]:
    """Read a gene-group file.

    Parameters
    ----------
    path : str | Path
        File with ``#`` metadata lines and the columns ``gene_id`` and ``symbol``.

    Returns
    -------
    tuple[pl.DataFrame, dict[str, str]]
        The genes and the metadata.

    Raises
    ------
    ValueError
        If the columns are missing, an ID is empty or duplicated.
    """
    metadata: dict[str, str] = {}
    rows: list[str] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            key, sep, value = line[1:].partition(":")
            if sep:
                metadata[key.strip()] = value.strip()
        elif line.strip():
            rows.append(line)
    if not rows or rows[0].split("\t")[:2] != ["gene_id", "symbol"]:
        raise ValueError(f"{path}: expected a header 'gene_id<TAB>symbol'")
    records = [dict(zip(("gene_id", "symbol"), r.split("\t")[:2], strict=True)) for r in rows[1:]]
    table = pl.DataFrame(records, schema={"gene_id": pl.Utf8, "symbol": pl.Utf8})
    if table.filter(pl.col("gene_id").str.strip_chars() == "").height:
        raise ValueError(f"{path}: empty gene_id")
    if table["gene_id"].n_unique() != table.height:
        raise ValueError(f"{path}: duplicated gene_id")
    return table, metadata


def default_gene_group(species: str, group: str) -> list[str]:
    """Return the stable Ensembl IDs of a shipped list.

    Raises
    ------
    FileNotFoundError
        If the list is not shipped (run ``build_gene_groups.py``).
    """
    path = gene_group_path(species, group)
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} not found; build it with recipes/common/build_gene_groups.py"
        )
    table, _ = read_gene_group(path)
    return table["gene_id"].to_list()


def excluded_symbols(symbols: list[str]) -> list[str]:
    """Return the symbols that ``SC-QC-006`` excludes from the groups (e.g. RPS6KA1, HBEGF)."""
    return [s for s in symbols if any(pattern.search(s) for pattern in _EXCLUDED)]
