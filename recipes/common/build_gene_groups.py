"""Build the default ribosomal and haemoglobin gene lists (``SC-QC-006``).

Usage (from the repository root)::

    python -m recipes.common.build_gene_groups --overwrite [--out DIR] [--snapshot-dir DIR]
    python -m recipes.common.build_gene_groups --from-snapshot DIR [--out DIR] [--overwrite]

Human: genes of the HGNC gene groups matched by ``HUMAN_GROUPS`` in the HGNC REST response for
protein-coding genes (approved, with an Ensembl ID), minus ``EXCLUDED_SYMBOL_PATTERNS``. Mouse:
protein-coding Ensembl genes of the MGI orthologues that HGNC records for those genes (BioMart).
Each file records the source URL, the download date, the SHA-256 of the HGNC response and its last
modification date (HGNC publishes no checksum), the Ensembl release and the SHA-256 of the source
snapshot. The builder fails when a group matches no gene or an excluded symbol remains, listing
candidate group names.

Sources: a build needs network access, and saves every response it uses, byte for byte, in a
snapshot directory (default ``build/gene_groups_sources/<date>``, which git ignores), with
``SOURCES.json`` (URL, parameters, size and SHA-256 of each file, download date) and
``SHA256SUMS``.
``--from-snapshot`` rebuilds the lists offline from such a directory after checking every checksum;
the lists are then identical to those of the original build. Keep the snapshot with the project's
inputs if the lists must be rebuildable. Outputs are written only after every source has been read,
each through a temporary file and ``os.replace``; existing lists (and a non-empty snapshot
directory) are replaced only with ``--overwrite`` (``CODE-PROV-005``).
"""

import argparse
import hashlib
import io
import json
import os
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
SPECIES = ("human", "mouse")
OUTPUT_NAMES = tuple(f"{species}_{group}.tsv" for group in HUMAN_GROUPS for species in SPECIES)
DEFAULT_SNAPSHOT_ROOT = Path("build") / "gene_groups_sources"
SNAPSHOT_INDEX = "SOURCES.json"
SNAPSHOT_SUMS = "SHA256SUMS"


def _get(url: str, **kwargs):
    import requests  # only the builder needs network access; not a package dependency

    response = requests.get(url, timeout=TIMEOUT, **kwargs)
    response.raise_for_status()
    return response


def write_atomic(path: Path, data: bytes) -> None:
    """Write ``data`` to ``path`` through a temporary file in its directory and ``os.replace``.

    A reader never sees a partly written file, and a failed write leaves an existing file intact.
    """
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        tmp.write_bytes(data)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sums_text(index: dict) -> str:
    """``SHA256SUMS`` content for the files of a snapshot index, sorted by file name."""
    entries = sorted(index["files"].values(), key=lambda entry: entry["file"])
    return "".join(f"{entry['sha256']}  {entry['file']}\n" for entry in entries)


def read_snapshot(directory: Path) -> dict:
    """Return the index of a source snapshot after checking every file against its SHA-256.

    Raises
    ------
    ValueError
        If ``SHA256SUMS`` disagrees with ``SOURCES.json``, a file name is not a plain name, or a
        file differs from its recorded checksum.
    """
    index = json.loads((directory / SNAPSHOT_INDEX).read_text(encoding="utf-8"))
    if (directory / SNAPSHOT_SUMS).read_text(encoding="utf-8") != _sums_text(index):
        raise ValueError(f"{directory}: {SNAPSHOT_SUMS} does not match {SNAPSHOT_INDEX}")
    for entry in index["files"].values():
        if Path(entry["file"]).name != entry["file"]:
            raise ValueError(f"{directory}: {entry['file']!r} is not a file name")
        if _sha256((directory / entry["file"]).read_bytes()) != entry["sha256"]:
            raise ValueError(f"{directory / entry['file']}: SHA-256 differs from {SNAPSHOT_INDEX}")
    return index


class Sources:
    """The responses of one build: fetched and saved in a snapshot, or read back from one.

    Parameters
    ----------
    directory
        Snapshot directory.
    offline
        Read the responses from ``directory`` (checked with ``read_snapshot``) instead of fetching.
    """

    def __init__(self, directory: Path, offline: bool = False):
        self.directory = directory
        self.offline = offline
        if offline:
            self.index = read_snapshot(directory)
        else:
            self.index = {"downloaded": date.today().isoformat(), "files": {}}

    @property
    def downloaded(self) -> str:
        """Date of the download (of the original build, for a rebuild from a snapshot)."""
        return self.index["downloaded"]

    def get(self, file: str, url: str, **kwargs) -> bytes:
        """Body of a GET request, saved as ``file`` in the snapshot (or read from it)."""
        params = kwargs.get("params")
        if self.offline:
            entry = self.index["files"].get(file)
            if entry is None or entry["url"] != url or entry["params"] != params:
                raise ValueError(f"{self.directory}: no saved response {file!r} for this request")
            return (self.directory / file).read_bytes()
        data = _get(url, **kwargs).content
        write_atomic(self.directory / file, data)
        self.index["files"][file] = {
            "file": file,
            "url": url,
            "params": params,
            "headers": kwargs.get("headers"),
            "bytes": len(data),
            "sha256": _sha256(data),
        }
        return data

    def sha256(self, file: str) -> str:
        """SHA-256 of one saved response."""
        return self.index["files"][file]["sha256"]

    def finish(self) -> str:
        """Write ``SOURCES.json`` and ``SHA256SUMS`` (except on a rebuild); return its SHA-256.

        The SHA-256 returned is that of ``SHA256SUMS``, recorded in every list.
        """
        sums = _sums_text(self.index)
        if not self.offline:
            index = json.dumps(self.index, indent=2, sort_keys=True) + "\n"
            write_atomic(self.directory / SNAPSHOT_INDEX, index.encode("utf-8"))
            write_atomic(self.directory / SNAPSHOT_SUMS, sums.encode("utf-8"))
        return _sha256(sums.encode("utf-8"))


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


def mouse_genes(
    human: pl.DataFrame, sources: Sources | None = None, file: str = "biomart.tsv"
) -> pl.DataFrame:
    """Protein-coding mouse genes of the MGI orthologues that HGNC records for the human genes.

    Ensembl Compara has no one-to-one orthologue for several core ribosomal proteins (e.g. RPL5,
    RPS6), because of their many retrocopies in mouse; the HGNC curation (``mgd_id``) has them.
    MGI IDs are mapped to Ensembl mouse genes with one BioMart query, saved as ``file`` in the
    snapshot of ``sources`` (fetched without saving when ``sources`` is None).
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
    params = {"query": query}
    if sources is None:
        body = _get(BIOMART, params=params).content
    else:
        body = sources.get(file, BIOMART, params=params)
    text = body.decode("utf-8")
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
    """Write one list with its metadata header, atomically (``write_atomic``)."""
    lines = [f"# {key}: {value}" for key, value in metadata.items()]
    lines.append("gene_id\tsymbol")
    lines += [f"{g}\t{s}" for g, s in table.select("gene_id", "symbol").iter_rows()]
    write_atomic(path, ("\n".join(lines) + "\n").encode("utf-8"))  # LF on every platform


def build(sources: Sources) -> dict[str, tuple[pl.DataFrame, dict[str, str]]]:
    """Read every source and return the lists by output file name, with their metadata."""
    headers = {"Accept": "application/json"}
    hgnc_file = "hgnc_protein_coding.json"
    hgnc = hgnc_table(json.loads(sources.get(hgnc_file, HGNC_URL, headers=headers)))
    info = json.loads(sources.get("hgnc_info.json", f"{HGNC_REST}/info", headers=headers))
    release = json.loads(
        sources.get(
            "ensembl_info_data.json",
            f"{ENSEMBL_REST}/info/data",
            headers={"Content-Type": "application/json"},
        )
    )
    tables = {}
    for group, table in human_groups(hgnc).items():
        tables[group] = (table, mouse_genes(table, sources, f"biomart_{group}.tsv"))
    sources_sha256 = sources.finish()
    lists = {}
    for group, (human, mouse) in tables.items():
        common = {
            "group": group,
            "source": HGNC_URL,
            "hgnc_groups": HUMAN_GROUPS[group],
            "hgnc_sha256": sources.sha256(hgnc_file),
            "hgnc_last_modified": info["lastModified"],
            "downloaded": sources.downloaded,
            "ensembl_release": str(release["releases"][0]),
            "sources_sha256": sources_sha256,
        }
        lists[f"human_{group}.tsv"] = (human, {"species": "human", **common})
        mouse_meta = {"species": "mouse", **common, "method": MOUSE_METHOD}
        lists[f"mouse_{group}.tsv"] = (mouse, mouse_meta)
    return lists


def main(argv: list[str] | None = None) -> int:
    """Read the sources (live or from a snapshot) and write the four lists."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=Path("recipes/common/data/gene_groups"))
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--snapshot-dir",
        type=Path,
        help=f"where to save the source responses (default {DEFAULT_SNAPSHOT_ROOT}/<date>)",
    )
    source.add_argument(
        "--from-snapshot",
        type=Path,
        help="rebuild offline from the responses saved in this directory",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="replace existing lists and a non-empty snapshot"
    )
    args = parser.parse_args(argv)
    existing = [str(args.out / name) for name in OUTPUT_NAMES if (args.out / name).exists()]
    if existing and not args.overwrite:
        print(
            f"ERROR: {', '.join(existing)} exist; pass --overwrite to replace them",
            file=sys.stderr,
        )
        return 1
    if args.from_snapshot:
        sources = Sources(args.from_snapshot, offline=True)
    else:
        snapshot = args.snapshot_dir or DEFAULT_SNAPSHOT_ROOT / date.today().isoformat()
        if snapshot.is_dir() and any(snapshot.iterdir()) and not args.overwrite:
            print(
                f"ERROR: {snapshot} is not empty; pass --overwrite to replace it", file=sys.stderr
            )
            return 1
        snapshot.mkdir(parents=True, exist_ok=True)
        sources = Sources(snapshot)
    lists = build(sources)
    args.out.mkdir(parents=True, exist_ok=True)
    for name, (table, metadata) in lists.items():
        write_group(args.out / name, table, metadata)
    for group in HUMAN_GROUPS:
        human, mouse = lists[f"human_{group}.tsv"][0], lists[f"mouse_{group}.tsv"][0]
        print(f"{group}: {human.height} human, {mouse.height} mouse genes")
    sums_sha256 = lists[OUTPUT_NAMES[0]][1]["sources_sha256"]
    print(f"sources: {sources.directory} ({SNAPSHOT_SUMS} SHA-256 {sums_sha256})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
