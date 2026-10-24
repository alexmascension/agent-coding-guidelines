"""Stable rule identifiers: assign them to new rules and keep the registry ``RULE_IDS.tsv``.

A rule is a top-level list item (``* `` or ``1. `` at column 0) of a guideline section, outside
code blocks and
tables. Its identifier is ``PREFIX-TOPIC-NNN`` (e.g. ``SC-QC-004``), written as inline code right
after the list
marker. Identifiers are never renumbered or reused: a new rule takes the next free number of its
topic, and a
removed rule's identifier is marked ``retired`` in the registry.

Usage (from the repository root)::

    python tools/rule_ids.py assign   # give identifiers to rules without one; update RULE_IDS.tsv
    python tools/rule_ids.py check    # verify identifiers against the guides and the registry
    (used in CI)
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "RULE_IDS.tsv"
GUIDES = ROOT / "guidelines"
ID_PATTERN = re.compile(r"`([A-Z]+-[A-Z]+-\d{3})`")
ITEM = re.compile(r"^(\* |\d+\. )(`[A-Z]+-[A-Z]+-\d{3}` )?")
SECTION = re.compile(r"^## (\d+)\. ")
SUBSECTION = re.compile(r"^### (.+)$")

PREFIX = {
    "CODING_GUIDELINES.md": "CODE",
    "STATISTICS_GUIDELINES.md": "STAT",
    "RNASEQ_BULK_GUIDELINES.md": "BULK",
    "SCRNASEQ_GUIDELINES.md": "SC",
    "SPATIAL_GUIDELINES.md": "SP",
    "META_ANALYSIS_GUIDELINES.md": "META",
    "SEX_INFERENCE_GUIDELINES.md": "SEX",
    "SAMPLE_IDENTITY_GUIDELINES.md": "ID",
    "DECONVOLUTION_GUIDELINES.md": "DECON",
    "HPC_GUIDELINES.md": "HPC",
}
# Topic code per section number. Sections absent from a map contain no rules: scope and audience
# notes (§0 of the
# domain guides), glossary, file lists, implementation tables.
TOPIC = {
    "CODE": {
        0: "USE",
        2: "PROJ",
        3: "LANG",
        4: "PRIN",
        5: "CHG",
        6: "NAME",
        7: "DATA",
        8: "TYPE",
        9: "DOC",
        10: "ROB",
        11: "NUM",
        12: "EXT",
        13: "PROV",
        14: "CONF",
        15: "REPRO",
        16: "TEST",
        17: "STYLE",
        18: "DEP",
        19: "PERF",
        20: "LOG",
        21: "IFACE",
        22: "NB",
        23: "VCS",
        24: "VAL",
        25: "AGENT",
    },
    "STAT": {
        0: "USE",
        1: "DES",
        2: "EXPL",
        3: "MOD",
        4: "HIER",
        5: "BATCH",
        6: "CIRC",
        7: "MT",
        8: "EFF",
        9: "ENR",
        10: "RES",
        11: "BAYES",
        12: "REP",
        13: "IMPL",
    },
    "BULK": {
        1: "DES",
        2: "REF",
        3: "RQC",
        4: "QUANT",
        5: "NORM",
        6: "EXPL",
        7: "HF",
        8: "DE",
        9: "TX",
        10: "ENR",
        11: "SPEC",
        12: "REP",
    },
    "SC": {
        1: "DATA",
        2: "DES",
        3: "QUANT",
        4: "QC",
        5: "NORM",
        6: "PC",
        7: "INT",
        8: "CLUST",
        9: "ANNOT",
        10: "DE",
        11: "DA",
        12: "CCC",
        13: "TRAJ",
        14: "SPEC",
        15: "REP",
    },
    "SP": {
        1: "PLAT",
        2: "DES",
        3: "DATA",
        4: "QC",
        5: "SEG",
        6: "NORM",
        7: "TYPE",
        8: "SVG",
        9: "DOM",
        10: "NBH",
        11: "SECT",
        12: "REP",
    },
    "META": {
        1: "PROT",
        2: "EXTR",
        3: "EFF",
        4: "MOD",
        5: "HET",
        6: "PUB",
        7: "SENS",
        8: "REP",
        9: "STRAT",
        10: "DATA",
        11: "SC",
        12: "OREP",
    },
    "SEX": {
        1: "WHY",
        2: "MARK",
        3: "COUNT",
        4: "REF",
        5: "THR",
        6: "DEC",
        7: "LOY",
        8: "SC",
        9: "MOUSE",
        10: "REP",
    },
    "ID": {1: "WHY", 2: "ORDER", 3: "DUP", 4: "GENO", 5: "MIX", 6: "POOL", 7: "RES", 8: "REP"},
    "DECON": {1: "EST", 2: "REFR", 3: "PREP", 4: "METH", 5: "VAL", 6: "DOWN", 7: "REP", 8: "SPAT"},
    "HPC": {
        1: "STOR",
        2: "FLOW",
        3: "ENV",
        4: "JOB",
        5: "WF",
        6: "IO",
        7: "PROV",
        8: "AGENT",
        9: "REP",
    },
}
# Subsections with their own topic (overrides the section topic from the subsection heading on).
SUBTOPIC = {
    ("CODE", "Recipes (reference implementations)"): "RCP",
    ("CODE", "Downloads from public archives"): "DL",
    ("HPC", "Datasets larger than the work area"): "BATCH",
}


def rule_lines(path: Path) -> list[tuple[int, str | None, str | None]]:
    """Return (line index, topic, existing ID) for every rule line of a guide."""
    prefix = PREFIX[path.name]
    topic = None
    in_code = False
    result = []
    for i, line in enumerate(path.read_text(encoding="utf-8").split("\n")):
        if line.startswith("```"):
            in_code = not in_code
            continue
        if in_code:
            continue
        section = SECTION.match(line)
        if section:
            topic = TOPIC[prefix].get(int(section.group(1)))
            continue
        if line.startswith("## "):
            topic = None  # part headings without number
            continue
        sub = SUBSECTION.match(line)
        if sub and (prefix, sub.group(1).strip()) in SUBTOPIC:
            topic = SUBTOPIC[(prefix, sub.group(1).strip())]
            continue
        item = ITEM.match(line)
        if item and topic is not None:
            existing = item.group(2).strip().strip("`") if item.group(2) else None
            result.append((i, topic, existing))
    return result


def read_registry() -> dict[str, tuple[str, str, str]]:
    """Return id -> (status, guide, summary)."""
    if not REGISTRY.exists():
        return {}
    rows = {}
    for line in REGISTRY.read_text(encoding="utf-8").splitlines()[1:]:
        rule_id, status, guide, summary = line.split("\t", 3)
        rows[rule_id] = (status, guide, summary)
    return rows


def summary_of(line: str) -> str:
    """Return the first words of a rule line, without markup, for the registry."""
    text = ITEM.sub("", line)
    text = re.sub(r"[*`]", "", text).strip()
    return text[:80]


def assign() -> None:
    """Give identifiers to rules without one; rewrite the registry, keeping retired IDs."""
    registry = read_registry()
    used = set(registry)
    rows = {k: v for k, v in registry.items() if v[0] == "retired"}
    for guide in PREFIX:
        path = GUIDES / guide
        lines = path.read_text(encoding="utf-8").split("\n")
        for i, topic, existing in rule_lines(path):
            if existing is None:
                stem = f"{PREFIX[guide]}-{topic}-"
                number = 1 + max((int(u[-3:]) for u in used if u.startswith(stem)), default=0)
                existing = f"{stem}{number:03d}"
                used.add(existing)
                marker = ITEM.match(lines[i]).group(1)
                lines[i] = f"{marker}`{existing}` {lines[i][len(marker) :]}"
            rows[existing] = ("active", guide, summary_of(lines[i]))
        path.write_text("\n".join(lines), encoding="utf-8", newline="\n")
    # identifiers in the old registry that are no longer in any guide become retired
    for rule_id, (_status, guide, summary) in registry.items():
        if rule_id not in rows:
            rows[rule_id] = ("retired", guide, summary)
    body = "\n".join(f"{k}\t{s}\t{g}\t{m}" for k, (s, g, m) in sorted(rows.items()))
    REGISTRY.write_text(
        "id\tstatus\tguide\tsummary\n" + body + "\n", encoding="utf-8", newline="\n"
    )


def check() -> list[str]:
    """Return problems: missing or duplicated IDs, wrong prefix, registry mismatch, reused IDs."""
    registry = read_registry()
    problems, seen = [], {}
    for guide in PREFIX:
        path = GUIDES / guide
        for i, _topic, existing in rule_lines(path):
            where = f"{guide}:{i + 1}"
            if existing is None:
                problems.append(f"{where}: rule without identifier (run tools/rule_ids.py assign)")
                continue
            if not existing.startswith(f"{PREFIX[guide]}-"):
                problems.append(f"{where}: {existing} has the prefix of another guide")
            if existing in seen:
                problems.append(f"{where}: {existing} duplicated (also {seen[existing]})")
            seen[existing] = where
            status = registry.get(existing, (None,))[0]
            if status is None:
                problems.append(f"{where}: {existing} not in RULE_IDS.tsv")
            elif status == "retired":
                problems.append(f"{where}: {existing} is retired and MUST NOT be reused")
    for rule_id, (status, guide, _) in registry.items():
        if status == "active" and rule_id not in seen:
            problems.append(
                f"RULE_IDS.tsv: {rule_id} active but absent from {guide} (mark it retired)"
            )
    return problems


def main() -> int:
    """Run ``assign`` or ``check``."""
    mode = sys.argv[1] if len(sys.argv) > 1 else "check"
    if mode == "assign":
        assign()
        return 0
    problems = check()
    for problem in problems:
        print(problem)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
