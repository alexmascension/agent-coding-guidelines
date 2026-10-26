r"""Check the guideline set: citation keys, cross-references, recipe paths, versions, BibTeX files.

Run from the repository root: ``python tools/check_guides.py``. Exits with status 1 and a list of
problems if any
check fails. Used in CI.

Checks
------
* Every ``[Author Year]`` key in an omics guide is a bold key of ``OMICS_REFERENCES.md``; every
  ``[Author Year†]``
  key, and every key in ``STATISTICS_GUIDELINES.md``, is a bold key of
  ``STATISTICS_REFERENCES.md``.
* Every ```X_GUIDELINES.md` §N`` reference points to an existing section of an existing guide.
* Guides (``guidelines/*_GUIDELINES.md``) and site profiles (``guidelines/hpc_sites/*.md``) are
  checked for both; reference files are in ``references/``.
* Every ``recipes/...`` path cited in a guide exists.
* Every guide declares ``Version MAJOR.MINOR.PATCH`` in its header.
* Every ``references/*_REFERENCES.bib`` has unique keys, no empty author names (Crossref exports a
  group author as an empty name) and months written as BibTeX macros (``sep``, not ``Sept``); it
  holds the same DOIs as its Markdown list, one entry per DOI on each side.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KEY = re.compile(r"[A-ZÀ-Ž][\w'’\-À-ž &]+ \d{4}[a-z]?†?")
GROUP = re.compile(r"\[([^\]\[]*?\d{4}[a-z]?†?[^\]\[]*?)\]")
XREF = re.compile(r"`(\w+_GUIDELINES\.md)` §(\d+(?:\.\d+)*)")
HEADING = re.compile(r"^#{2,4} (\d+(?:\.\d+)*)\.? ", re.M)
RECIPE = re.compile(r"`(recipes/[\w/]+\.(?:py|R|sh|slurm|env))`")
VERSION = re.compile(r"Version \d+\.\d+\.\d+")
BOLD_KEY = re.compile(r"^- \*\*\[([^\]]+)\]\*\*", re.M)
BIB_KEY = re.compile(r"^@\w+\{([^,\s]+),", re.M)
BIB_DOI = re.compile(r"\bDOI=\{([^}]+)\}", re.I)
MD_DOI = re.compile(r"^- \*\*\[[^\]]+\]\*\*.*?\[doi:([^\]]+)\]", re.M)
EMPTY_AUTHOR = re.compile(r"\bauthor=\{\s*and\s|\sand\s+and\s|\sand\s*\}")
BIB_MONTH = re.compile(r"\bmonth=([A-Za-z]+)\b")
MONTH_MACROS = {"jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"}


def bold_keys(path: Path) -> set[str]:
    """Return the bold keys defined in a references file."""
    return set(BOLD_KEY.findall(path.read_text(encoding="utf-8"))) if path.exists() else set()


def repeated(items: list[str]) -> list[str]:
    """Return the items that occur more than once, each once, in order of first repetition."""
    seen, found = set(), []
    for item in items:
        if item in seen and item not in found:
            found.append(item)
        seen.add(item)
    return found


def bib_problems(name: str, bib: str, md: str) -> list[str]:
    """Return the problems of a BibTeX file and of its agreement with its Markdown list."""
    keys = BIB_KEY.findall(bib)
    problems = [f"{name}.bib: key {key} used more than once" for key in repeated(keys)]
    bib_dois = [doi.lower() for doi in BIB_DOI.findall(bib)]
    md_dois = [doi.lower() for doi in MD_DOI.findall(md)]
    problems += [f"{name}.bib: DOI {doi} in more than one entry" for doi in repeated(bib_dois)]
    problems += [f"{name}.md: DOI {doi} in more than one entry" for doi in repeated(md_dois)]
    for doi in sorted(set(md_dois) - set(bib_dois)):
        problems.append(f"{name}.md: DOI {doi} has no entry in {name}.bib")
    for doi in sorted(set(bib_dois) - set(md_dois)):
        problems.append(f"{name}.bib: DOI {doi} has no entry in {name}.md")
    for number, line in enumerate(bib.splitlines(), start=1):
        if EMPTY_AUTHOR.search(line):
            problems.append(f"{name}.bib:{number}: empty author name")
        for month in BIB_MONTH.findall(line):
            if month.lower() not in MONTH_MACROS:
                problems.append(f"{name}.bib:{number}: month={month} is not a BibTeX month macro")
    return problems


def main() -> int:
    """Run all checks and print the problems found."""
    guides = sorted((ROOT / "guidelines").glob("*_GUIDELINES.md"))
    omics_keys = bold_keys(ROOT / "references" / "OMICS_REFERENCES.md")
    stats_keys = bold_keys(ROOT / "references" / "STATISTICS_REFERENCES.md")
    headings = {g.name: set(HEADING.findall(g.read_text(encoding="utf-8"))) for g in guides}
    problems = []
    for guide in [*guides, *sorted((ROOT / "guidelines" / "hpc_sites").glob("*.md"))]:
        text = guide.read_text(encoding="utf-8")
        if guide in guides and not VERSION.search(text[:500]):
            problems.append(f"{guide.name}: no version in header")
        for group in GROUP.findall(text):
            for key in (k.strip() for k in group.split(";")):
                if not KEY.fullmatch(key) or guide.name == "CODING_GUIDELINES.md":
                    continue
                if key.endswith("†") or guide.name == "STATISTICS_GUIDELINES.md":
                    if key.rstrip("†") not in stats_keys:
                        problems.append(f"{guide.name}: [{key}] not in STATISTICS_REFERENCES.md")
                elif key not in omics_keys:
                    problems.append(f"{guide.name}: [{key}] not in OMICS_REFERENCES.md")
        for target, section in XREF.findall(text):
            if target not in headings:
                problems.append(f"{guide.name}: reference to missing guide {target}")
            elif section not in headings[target]:
                problems.append(f"{guide.name}: {target} §{section} does not exist")
        for path in RECIPE.findall(text):
            if not (ROOT / path).exists():
                problems.append(f"{guide.name}: recipe {path} does not exist")
    for bib_path in sorted((ROOT / "references").glob("*_REFERENCES.bib")):
        md_path = bib_path.with_suffix(".md")
        md = md_path.read_text(encoding="utf-8") if md_path.exists() else ""
        problems += bib_problems(bib_path.stem, bib_path.read_text(encoding="utf-8"), md)
    for problem in problems:
        print(problem)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
