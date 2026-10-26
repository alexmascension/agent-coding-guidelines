"""Decide which CI test suites a change needs: Python recipe tests and R-backed tests.

Run from the repository root::

    python tools/ci_scope.py --base <ref> [--head <ref>] [--github]

Prints ``python=true|false`` and ``r=true|false`` (and appends them to ``$GITHUB_OUTPUT`` with
``--github``). Without a usable base (empty, all zeros, or unknown to git: scheduled or manual
runs, new branches) every suite runs.

Rules
-----
* ``python``: a file under ``recipes/`` or ``tools/`` changed, or ``pyproject.toml``, or the CI
  definition.
* ``r``: an ``*.R`` file changed; or a changed file depends on R at the base or at the head, so a
  change that adds, removes or alters a dependency on R is caught. A file depends on R when it
  calls R or names an R script (``Rscript``, ``RECIPES_RSCRIPT``, ``<name>.R``), when it is a
  module used by a test marked ``r`` (e.g. the Python implementation compared with R), or when
  such a file imports it, directly or transitively; or ``pyproject.toml`` or the CI definition
  changed.
"""

import argparse
import ast
import os
import re
import subprocess
import sys
from collections.abc import Callable, Iterable
from pathlib import PurePosixPath

CI_FILES = {".github/workflows/ci.yml", "tools/ci_scope.py"}
CALLS_R = re.compile(r"\bRscript\b|RECIPES_RSCRIPT|[\w/]+\.R\b")
EMPTY = {"", "0" * 40}


def git(*args: str) -> str:
    """Run git and return its standard output."""
    return subprocess.run(["git", *args], capture_output=True, text=True, check=True).stdout


def module_path(name: str, files: set[str]) -> str | None:
    """Return the file of a ``recipes`` module name, or None when it is not in the tree."""
    stem = name.replace(".", "/")
    for path in (f"{stem}.py", f"{stem}/__init__.py"):
        if path in files:
            return path
    return None


def imported_modules(path: str, text: str, files: set[str]) -> set[str]:
    """Return the files of the ``recipes`` modules imported by a Python file."""
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return set()
    package = str(PurePosixPath(path).parent).replace("/", ".")
    found = set()
    for node in ast.walk(tree):
        names = []
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            base = node.module or ""
            if node.level:
                parts = package.split(".")[: len(package.split(".")) - node.level + 1]
                base = ".".join([*parts, base] if base else parts)
            names = [base, *(f"{base}.{alias.name}" for alias in node.names)]
        for name in names:
            if name.startswith("recipes") and (target := module_path(name, files)):
                found.add(target)
    return found


def r_test_modules(path: str, text: str, files: set[str]) -> set[str]:
    """Return the ``recipes`` modules used by the tests marked ``r`` in a test file."""
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return set()
    bound = {}
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("recipes"):
            for alias in node.names:
                bound[alias.asname or alias.name] = node.module
    used = set()
    for node in tree.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        marked = any(
            "mark.r" == ast.unparse(d).removeprefix("pytest.") for d in node.decorator_list
        )
        if not marked:
            continue
        used |= imported_modules(path, ast.unparse(node), files)
        for name in ast.walk(node):
            if isinstance(name, ast.Name) and name.id in bound:
                if target := module_path(bound[name.id], files):
                    used.add(target)
    return used


def r_dependent(files: Iterable[str], read: Callable[[str], str | None]) -> set[str]:
    """Return the files that depend on R or that R-backed code and tests depend on.

    Seeds: modules and scripts that call R or name an R script, test files with tests marked
    ``r``, and the modules those tests use. The result adds the ``recipes`` modules the seeds
    import, transitively. Other imports of a test file do not count.
    """
    files = set(files)
    texts = {f: read(f) for f in files if f.endswith((".py", ".sh", ".slurm", ".toml"))}
    seeds = set()
    for path, text in texts.items():
        if not text:
            continue
        if "/tests/" in path:
            modules = r_test_modules(path, text, files)
            if modules or CALLS_R.search(text):
                seeds |= {path, *modules}
        elif CALLS_R.search(text):
            seeds.add(path)
    closure, queue = set(seeds), [s for s in seeds if "/tests/" not in s]
    while queue:
        current = queue.pop()
        if current.endswith(".py") and texts.get(current):
            for target in imported_modules(current, texts[current], files) - closure:
                closure.add(target)
                queue.append(target)
    return closure


def scope(changed: set[str], trees: list[tuple[set[str], Callable[[str], str | None]]]) -> dict:
    """Return which suites the changed files need, given the base and head trees."""
    shared = bool(changed & (CI_FILES | {"pyproject.toml"}))
    python = shared or any(f.startswith(("recipes/", "tools/")) for f in changed)
    r = shared or any(f.endswith(".R") for f in changed)
    for files, read in trees:
        r = r or bool(changed & r_dependent(files, read))
    return {"python": python, "r": r}


def tree_at(ref: str) -> tuple[set[str], Callable[[str], str | None]]:
    """Return the files under ``recipes/`` at a git ref and a reader for their contents."""
    files = set(git("ls-tree", "-r", "--name-only", ref, "--", "recipes").split())

    def read(path: str) -> str | None:
        result = subprocess.run(["git", "show", f"{ref}:{path}"], capture_output=True, text=True)
        return result.stdout if result.returncode == 0 else None

    return files, read


def main() -> int:
    """Print the suites to run and, with ``--github``, write them as step outputs."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--base", default="")
    parser.add_argument("--head", default="HEAD")
    parser.add_argument("--github", action="store_true")
    args = parser.parse_args()
    known = args.base not in EMPTY and (
        subprocess.run(["git", "cat-file", "-e", f"{args.base}^{{commit}}"]).returncode == 0
    )
    if known:
        changed = set(git("diff", "--name-only", args.base, args.head).split())
        result = scope(changed, [tree_at(args.base), tree_at(args.head)])
        print(f"changed files: {len(changed)}", file=sys.stderr)
    else:
        result = {"python": True, "r": True}
        print("no usable base: running every suite", file=sys.stderr)
    lines = [f"{key}={str(value).lower()}" for key, value in result.items()]
    print("\n".join(lines))
    if args.github and os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
