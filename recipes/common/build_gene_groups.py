"""Build the default ribosomal and haemoglobin gene lists (``SC-QC-006``); needs network access.

Usage (from the repository root)::

    python -m recipes.common.build_gene_groups [--out recipes/common/data/gene_groups]

Human: genes of the HGNC gene groups matched by ``HUMAN_GROUPS`` in the HGNC REST response for
protein-coding genes (approved, with an Ensembl ID), minus ``EXCLUDED_SYMBOL_PATTERNS``. Mouse:
protein-coding Ensembl genes of the MGI orthologues that HGNC records for those genes (BioMart).
Each file records the source URL, the download date, the SHA-256 of the HGNC response and its last
modification date (HGNC publishes no checksum) and the Ensembl release. The builder fails when a
group matches no gene or an excluded symbol remains, listing candidate group names.
"""

import argparse
import hashlib
import io
import json
import re
import sys
from datetime import date
from pathlib import Path

import polars as pl

from recipes.common.gene_groups import excluded_symbols

HGNC_REST = "https://rest.genenames.org"
HGNC_URL = f"{HGNC_REST}/fetch/locus_group/protein-coding%20gene"
ENSEMBL_REST = "https://rest.ensembl.org"
BIOMART = "https://www.ensembl.org/biomart/martservice"
MOUSE_METHOD = "MGI orthologues curated by HGNC (mgd_id), Ensembl protein-coding genes via BioMart"
HUMAN_GROUPS = {
    "ribosomal": r"^(L|S) ribosomal proteins$",
    "haemoglobin": r"^Hemoglobin subunits$",
}
CANDIDATE_HINTS = {"ribosomal": "ribosom", "haemoglobin": "globin"}
TIMEOUT = 300  # seconds; the HGNC response for all protein-coding genes is about 19 MB


def _get(url: str, **kwargs):
    import requests  # only the builder needs network access; not a package dependency

    response = requests.get(url, timeout=TIMEOUT, **kwargs)
    response.raise_for_status()
    return response


def hgnc_table(response: dict) -> pl.DataFrame:
    """Table of the HGNC REST response, with gene groups joined by ``|`` as in the HGNC files."""
    rows = [
        {
            "symbol": doc.get("symbol", ""),
            "status": doc.get("status", ""),
            "locus_group": doc.get("locus_group", ""),
            "ensembl_gene_id": doc.get("ensembl_gene_id") or "",
            "gene_group": "|".join(doc.get("gene_group", [])),
            "mgd_id": "|".join(doc.get("mgd_id", [])),
        }
        for doc in response["response"]["docs"]
    ]
    return pl.DataFrame(rows, schema=dict.fromkeys(rows[0], pl.Utf8))


def human_groups(hgnc: pl.DataFrame) -> dict[str, pl.DataFrame]:
    """Select the genes of each group of ``HUMAN_GROUPS`` from the HGNC complete set."""
    approved = hgnc.filter(
        (pl.col("status") == "Approved")
        & (pl.col("locus_group") == "protein-coding gene")
        & pl.col("ensembl_gene_id").is_not_null()
        & (pl.col("ensembl_gene_id") != "")
    )
    out = {}
    for group, pattern in HUMAN_GROUPS.items():
        regex = re.compile(pattern)
        keep = [
            any(regex.match(name.strip()) for name in (groups or "").split("|"))
            for groups in approved["gene_group"].to_list()
        ]
        table = approved.filter(pl.Series(keep)).select(
            pl.col("ensembl_gene_id").alias("gene_id"), "symbol", "mgd_id"
        )
        table = table.filter(~pl.col("symbol").is_in(excluded_symbols(table["symbol"].to_list())))
        if table.height == 0:
            names = sorted(
                {
                    name.strip()
                    for groups in hgnc["gene_group"].drop_nulls().to_list()
                    for name in groups.split("|")
                    if CANDIDATE_HINTS[group] in name.lower()
                }
            )
            raise ValueError(f"no HGNC gene matches {pattern!r}; candidate groups: {names[:30]}")
        out[group] = table.sort("symbol")
    return out


def mouse_genes(human: pl.DataFrame) -> pl.DataFrame:
    """Protein-coding mouse genes of the MGI orthologues that HGNC records for the human genes.

    Ensembl Compara has no one-to-one orthologue for several core ribosomal proteins (e.g. RPL5,
    RPS6), because of their many retrocopies in mouse; the HGNC curation (``mgd_id``) has them.
    MGI IDs are mapped to Ensembl mouse genes with one BioMart query.
    """
    mgi = sorted({m for ids in human["mgd_id"].to_list() for m in ids.split("|") if m})
    if not mgi:
        raise ValueError("no MGI orthologue recorded by HGNC")
    query = (
        '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" '
        'formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">'
        '<Dataset name="mmusculus_gene_ensembl" interface="default">'
        f'<Filter name="mgi_id" value="{",".join(mgi)}"/>'
        '<Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/>'
        '<Attribute name="gene_biotype"/></Dataset></Query>'
    )
    text = _get(BIOMART, params={"query": query}).text
    table = pl.read_csv(io.StringIO(text), separator="\t", infer_schema_length=0)
    if table.columns != ["Gene stable ID", "Gene name", "Gene type"]:
        raise ValueError(f"unexpected BioMart response: {text[:200]!r}")
    table = (
        table.filter(pl.col("Gene type") == "protein_coding")
        .select(pl.col("Gene stable ID").alias("gene_id"), pl.col("Gene name").alias("symbol"))
        .unique()
        .sort("symbol")
    )
    return table.filter(~pl.col("symbol").is_in(excluded_symbols(table["symbol"].to_list())))


def write_group(path: Path, table: pl.DataFrame, metadata: dict[str, str]) -> None:
    """Write one list with its metadata header."""
    lines = [f"# {key}: {value}" for key, value in metadata.items()]
    lines.append("gene_id\tsymbol")
    lines += [f"{g}\t{s}" for g, s in table.select("gene_id", "symbol").iter_rows()]
    path.write_bytes(("\n".join(lines) + "\n").encode("utf-8"))  # LF on every platform


def main() -> int:
    """Download the sources and write the four lists."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="recipes/common/data/gene_groups")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    headers = {"Accept": "application/json"}
    raw = _get(HGNC_URL, headers=headers).content
    sha256 = hashlib.sha256(raw).hexdigest()
    hgnc = hgnc_table(json.loads(raw))
    hgnc_modified = _get(f"{HGNC_REST}/info", headers=headers).json()["lastModified"]
    release = _get(f"{ENSEMBL_REST}/info/data", headers={"Content-Type": "application/json"})
    ensembl_release = str(release.json()["releases"][0])
    today = date.today().isoformat()

    for group, table in human_groups(hgnc).items():
        common = {
            "group": group,
            "source": HGNC_URL,
            "hgnc_groups": HUMAN_GROUPS[group],
            "hgnc_sha256": sha256,
            "hgnc_last_modified": hgnc_modified,
            "downloaded": today,
            "ensembl_release": ensembl_release,
        }
        write_group(out / f"human_{group}.tsv", table, {"species": "human", **common})
        mouse = mouse_genes(table)
        mouse_meta = {
            "species": "mouse",
            **common,
            "method": MOUSE_METHOD,
        }
        write_group(out / f"mouse_{group}.tsv", mouse, mouse_meta)
        print(f"{group}: {table.height} human, {mouse.height} mouse genes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
