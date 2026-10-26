"""Report changes in the upstream repositories the guides adapt, since the commits last reviewed.

Run from the repository root::

    python tools/check_sources.py [--output report.md] [--github-output]
    python tools/check_sources.py mark-reviewed <name> <commit>

The registry ``references/SOURCE_REPOS.toml`` lists each upstream repository with the commit last
reviewed and the globs of the files followed in it. For each repository, the GitHub compare API
(``/repos/{owner}/{repo}/compare/{reviewed_commit}...{branch}``) gives the files changed since
that commit; each changed file takes the status of the first glob it matches. ``GITHUB_TOKEN``,
when set, authenticates the requests. A review is needed when a file with status ``integrated``,
``deferred`` or ``watch`` changed; ``--github-output`` writes ``changed=true|false`` for the
monthly workflow (``.github/workflows/sources.yml``). ``mark-reviewed`` records the commit a
review covered, in the block of that repository only.
"""

import argparse
import datetime
import fnmatch
import json
import os
import re
import sys
import tomllib
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "references" / "SOURCE_REPOS.toml"
API = "https://api.github.com"
COMPARE_FILE_LIMIT = 300  # the compare API lists at most 300 changed files
STATUSES = ("integrated", "deferred", "watch", "not-integrated")
REVIEW_ORDER = ("integrated", "deferred", "watch")  # a change to these opens a review
REVIEW_STATUSES = set(REVIEW_ORDER)
REPO_KEYS = {"name", "url", "branch", "licence", "reviewed_commit", "reviewed_on", "path"}
GITHUB_URL = re.compile(r"^https://github\.com/([\w.-]+/[\w.-]+)$")
COMMIT = re.compile(r"^[0-9a-f]{40}$")


@dataclass
class RepoReport:
    """Changes in one upstream repository since its reviewed commit."""

    name: str
    compare_url: str
    commits_ahead: int
    tracked: dict[str, list[tuple[str, str, str]]] = field(default_factory=dict)
    untracked: list[str] = field(default_factory=list)
    truncated: bool = False

    def needs_review(self) -> bool:
        """Return True when a file whose status calls for a review changed."""
        return any(self.tracked.get(status) for status in REVIEW_STATUSES)


def load_registry(path: Path = REGISTRY) -> list[dict]:
    """Return the repositories of the registry, after checking their fields and statuses."""
    with path.open("rb") as handle:
        repos = tomllib.load(handle).get("repo", [])
    for repo in repos:
        missing = REPO_KEYS - set(repo)
        if missing:
            raise ValueError(f"{repo.get('name', '?')}: missing keys {sorted(missing)}")
        owner_repo(repo["url"])
        if not COMMIT.match(repo["reviewed_commit"]):
            raise ValueError(f"{repo['name']}: reviewed_commit is not a full SHA")
        for entry in repo["path"]:
            if entry.get("status") not in STATUSES:
                raise ValueError(
                    f"{repo['name']}: {entry.get('glob')!r} has status {entry.get('status')!r}"
                )
    return repos


def owner_repo(url: str) -> str:
    """Return ``owner/repo`` from a GitHub repository URL."""
    match = GITHUB_URL.match(url)
    if not match:
        raise ValueError(f"Not a GitHub repository URL: {url!r}")
    return match.group(1)


def fetch_compare(repo: dict, token: str | None) -> dict:
    """Return the compare API response from the reviewed commit to the head of the branch."""
    span = f"{repo['reviewed_commit']}...{repo['branch']}"
    url = f"{API}/repos/{owner_repo(repo['url'])}/compare/{span}"
    headers = {"Accept": "application/vnd.github+json", "User-Agent": "agent-coding-guidelines"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    with urllib.request.urlopen(
        urllib.request.Request(url, headers=headers), timeout=60
    ) as response:
        return json.load(response)


def matching_path(filename: str, paths: list[dict]) -> dict | None:
    """Return the first followed path whose glob matches ``filename``."""
    for entry in paths:
        if fnmatch.fnmatchcase(filename, entry["glob"]):
            return entry
    return None


def summarize(repo: dict, compare: dict) -> RepoReport:
    """Group the files changed in a compare response by the status of the glob they match.

    A renamed file is matched by its new name, then by its previous one. Files that match no glob
    are listed as untracked.
    """
    report = RepoReport(
        name=repo["name"],
        compare_url=compare.get("html_url", repo["url"]),
        commits_ahead=compare.get("ahead_by", 0),
        tracked={status: [] for status in STATUSES},
    )
    files = compare.get("files", [])
    report.truncated = len(files) >= COMPARE_FILE_LIMIT
    for changed in files:
        entry = matching_path(changed["filename"], repo["path"])
        if entry is None and changed.get("previous_filename"):
            entry = matching_path(changed["previous_filename"], repo["path"])
        if entry is None:
            report.untracked.append(changed["filename"])
        else:
            change = (changed["filename"], changed["status"], entry.get("into", ""))
            report.tracked[entry["status"]].append(change)
    return report


def render(reports: list[RepoReport], today: datetime.date) -> str:
    """Return the Markdown report of all repositories."""
    lines = [f"# Upstream source review ({today.isoformat()})", ""]
    lines.append(
        "Changes since the commits recorded in `references/SOURCE_REPOS.toml`. Review each file "
        "listed under integrated, deferred or watch, adapt what is worth adapting, then run "
        "`python tools/check_sources.py mark-reviewed <name> <commit>` in the same pull request."
    )
    for report in reports:
        lines += ["", f"## {report.name}", ""]
        if report.commits_ahead == 0:
            lines.append("No new commits.")
            continue
        lines.append(f"{report.commits_ahead} new commits: <{report.compare_url}>")
        if report.truncated:
            lines.append(
                f"The compare API listed only {COMPARE_FILE_LIMIT} files; check the rest "
                "on the compare page."
            )
        for status in REVIEW_ORDER:
            for filename, change, into in report.tracked[status]:
                lines.append(f"- [{status}] `{filename}` ({change}): {into}")
        skipped = report.tracked["not-integrated"]
        if skipped:
            lines += ["", f"<details><summary>{len(skipped)} not-integrated files</summary>", ""]
            lines += [f"- `{filename}` ({change}): {into}" for filename, change, into in skipped]
            lines += ["", "</details>", ""]
        if report.untracked:
            lines.append(f"- {len(report.untracked)} other files changed, not followed.")
    return "\n".join(lines) + "\n"


def mark_reviewed(text: str, name: str, commit: str, day: datetime.date) -> str:
    """Return the registry text with the reviewed commit and date of one repository replaced."""
    if not COMMIT.match(commit):
        raise ValueError(f"Not a full commit SHA: {commit!r}")
    blocks = re.split(r"(?m)^(?=\[\[repo\]\]\s*$)", text)
    found = 0
    for index, block in enumerate(blocks):
        if not re.search(rf'(?m)^name = "{re.escape(name)}"\s*$', block):
            continue
        found += 1
        block, commits = re.subn(
            r'(?m)^reviewed_commit = "[^"]*"', f'reviewed_commit = "{commit}"', block
        )
        block, days = re.subn(r'(?m)^reviewed_on = "[^"]*"', f'reviewed_on = "{day}"', block)
        if commits != 1 or days != 1:
            raise ValueError(f"{name}: expected one reviewed_commit and one reviewed_on line")
        blocks[index] = block
    if found != 1:
        raise ValueError(f"{name!r} names {found} repositories in the registry, expected 1")
    return "".join(blocks)


def main() -> int:
    """Write the report, or record a reviewed commit."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command")
    marker = sub.add_parser("mark-reviewed", help="record the commit a review covered")
    marker.add_argument("name")
    marker.add_argument("commit")
    parser.add_argument("--output", type=Path, help="write the Markdown report to this file")
    parser.add_argument(
        "--github-output", action="store_true", help="append changed=true|false to $GITHUB_OUTPUT"
    )
    args = parser.parse_args()
    today = datetime.datetime.now(datetime.UTC).date()
    if args.command == "mark-reviewed":
        text = REGISTRY.read_text(encoding="utf-8")
        new_text = mark_reviewed(text, args.name, args.commit, today)
        REGISTRY.write_text(new_text, encoding="utf-8", newline="\n")
        return 0
    token = os.environ.get("GITHUB_TOKEN")
    reports = []
    for repo in load_registry():
        try:
            compare = fetch_compare(repo, token)
        except urllib.error.HTTPError as error:
            print(f"{repo['name']}: compare request failed: HTTP {error.code}", file=sys.stderr)
            return 1
        reports.append(summarize(repo, compare))
    report_text = render(reports, today)
    if args.output:
        args.output.write_text(report_text, encoding="utf-8", newline="\n")
    else:
        print(report_text, end="")
    changed = any(report.needs_review() for report in reports)
    if args.github_output:
        with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as handle:
            handle.write(f"changed={'true' if changed else 'false'}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
