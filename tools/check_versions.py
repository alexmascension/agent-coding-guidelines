"""Check that every changed guide bumps its version and has a changelog entry (CODE-USE-007).

Run from the repository root::

    python tools/check_versions.py [--base <ref>]

Compares the working tree with ``--base`` (default ``HEAD``; in CI, the commit before the push or
the base of the pull request). For each ``guidelines/*_GUIDELINES.md`` whose content differs from
the base (a file moved without changes is not a change), the
``Version MAJOR.MINOR.PATCH`` of its header MUST be greater than at the base, and ``CHANGELOG.md``
MUST contain the heading ``### <file> [<version>]``. Without a usable base the check is skipped.
Exits with status 1 and a list of problems if any check fails. Used in pre-commit and CI.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION = re.compile(r"Version (\d+)\.(\d+)\.(\d+)")
EMPTY = {"", "0" * 40}


def version(text: str) -> tuple[int, int, int] | None:
    """Return the version declared in the first lines of a guide."""
    match = VERSION.search(text[:500])
    return tuple(int(part) for part in match.groups()) if match else None


def main() -> int:
    """Run the check and print the problems found."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--base", default="HEAD")
    args = parser.parse_args()
    if (
        args.base in EMPTY
        or subprocess.run(
            ["git", "cat-file", "-e", f"{args.base}^{{commit}}"], cwd=ROOT
        ).returncode
    ):
        print("check_versions: no usable base, skipped")
        return 0
    status = subprocess.run(
        ["git", "diff", "--name-status", "-M", args.base, "--", "*_GUIDELINES.md"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.splitlines()
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    problems = []
    for line in status:
        code, *paths = line.split("\t")
        if code.startswith("D") or code == "R100":  # deleted, or moved without changes
            continue
        old_name, name = paths[0], paths[-1]
        path = ROOT / name
        if not path.exists():
            continue
        new = version(path.read_text(encoding="utf-8"))
        old_text = subprocess.run(
            ["git", "show", f"{args.base}:{old_name}"], cwd=ROOT, capture_output=True, text=True
        )
        old = version(old_text.stdout) if old_text.returncode == 0 else None
        shown = ".".join(map(str, new)) if new else "?"
        heading = f"### {Path(name).name} [{shown}]"
        if new is None:
            problems.append(f"{name}: no version in header")
        elif old is not None and new <= old:
            problems.append(f"{name}: changed but version {shown} is not greater than at base")
        if new is not None and heading not in changelog:
            problems.append(f"{name}: CHANGELOG.md has no heading '{heading}'")
    for problem in problems:
        print(problem)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
