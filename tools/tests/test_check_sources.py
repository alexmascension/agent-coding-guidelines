"""Tests of ``tools/check_sources.py`` on registry text and compare responses written by hand."""

import datetime

import pytest

from tools.check_sources import (
    COMPARE_FILE_LIMIT,
    load_registry,
    mark_reviewed,
    owner_repo,
    render,
    summarize,
)

SHA_A = "a" * 40
SHA_B = "b" * 40
SHA_C = "c" * 40
DAY = datetime.date(2026, 11, 1)

REGISTRY = f"""
[[repo]]
name = "alpha"
url = "https://github.com/owner/alpha"
branch = "main"
licence = "MIT"
reviewed_commit = "{SHA_A}"
reviewed_on = "2026-10-26"

[[repo.path]]
glob = "skills/brainstorming/SKILL.md"
status = "integrated"
into = "AGENTIC §13"

[[repo.path]]
glob = "skills/brainstorming/*"
status = "not-integrated"
into = "tooling"

[[repo.path]]
glob = "skills/tdd/*"
status = "deferred"
into = "CODING round"

[[repo.path]]
glob = "skills/*/SKILL.md"
status = "watch"
into = "new skills"

[[repo]]
name = "beta"
url = "https://github.com/owner/beta"
branch = "main"
licence = "MIT"
reviewed_commit = "{SHA_B}"
reviewed_on = "2026-10-26"

[[repo.path]]
glob = "guide.md"
status = "integrated"
into = "AGENTIC §1"
"""


@pytest.fixture
def repos(tmp_path):
    path = tmp_path / "SOURCE_REPOS.toml"
    path.write_text(REGISTRY, encoding="utf-8")
    return load_registry(path)


def compare(*files, ahead_by=1):
    return {
        "html_url": "https://github.com/owner/alpha/compare/aaa...main",
        "ahead_by": ahead_by,
        "files": [dict(item) for item in files],
    }


def test_groups_changed_files_by_the_first_matching_glob(repos):
    report = summarize(
        repos[0],
        compare(
            {"filename": "skills/brainstorming/SKILL.md", "status": "modified"},
            {"filename": "skills/brainstorming/server.cjs", "status": "modified"},
            {"filename": "skills/tdd/SKILL.md", "status": "modified"},
            {"filename": "skills/new-skill/SKILL.md", "status": "added"},
            {"filename": "tests/run.sh", "status": "modified"},
        ),
    )
    assert report.tracked["integrated"] == [
        ("skills/brainstorming/SKILL.md", "modified", "AGENTIC §13")
    ]
    assert report.tracked["not-integrated"] == [
        ("skills/brainstorming/server.cjs", "modified", "tooling")
    ]
    assert report.tracked["deferred"] == [("skills/tdd/SKILL.md", "modified", "CODING round")]
    assert report.tracked["watch"] == [("skills/new-skill/SKILL.md", "added", "new skills")]
    assert report.untracked == ["tests/run.sh"]
    assert report.needs_review()
    text = render([report], DAY)
    assert "- [integrated] `skills/brainstorming/SKILL.md` (modified): AGENTIC §13" in text
    assert "<details><summary>1 not-integrated files</summary>" in text
    assert "- 1 other files changed, not followed." in text


def test_matches_a_renamed_file_by_its_previous_name(repos):
    renamed = {
        "filename": "docs/brainstorming.md",
        "previous_filename": "skills/brainstorming/SKILL.md",
        "status": "renamed",
    }
    report = summarize(repos[0], compare(renamed))
    assert report.tracked["integrated"] == [("docs/brainstorming.md", "renamed", "AGENTIC §13")]


def test_changes_only_to_untracked_or_not_integrated_files_need_no_review(repos):
    report = summarize(
        repos[0],
        compare(
            {"filename": "tests/run.sh", "status": "modified"},
            {"filename": "skills/brainstorming/server.cjs", "status": "modified"},
        ),
    )
    assert not report.needs_review()


def test_no_new_commits_reports_nothing(repos):
    report = summarize(repos[1], compare(ahead_by=0))
    assert not report.needs_review()
    assert "No new commits." in render([report], DAY)


def test_flags_a_file_list_truncated_by_the_api(repos):
    files = [{"filename": f"other/{i}.txt", "status": "added"} for i in range(COMPARE_FILE_LIMIT)]
    report = summarize(repos[1], compare(*files, ahead_by=400))
    assert report.truncated
    assert f"listed only {COMPARE_FILE_LIMIT} files" in render([report], DAY)


def test_mark_reviewed_changes_only_the_named_repository():
    text = mark_reviewed(REGISTRY, "beta", SHA_C, DAY)
    assert f'reviewed_commit = "{SHA_A}"' in text
    assert f'reviewed_commit = "{SHA_C}"' in text
    assert f'reviewed_commit = "{SHA_B}"' not in text
    assert text.count('reviewed_on = "2026-11-01"') == 1
    assert text.count('reviewed_on = "2026-10-26"') == 1


def test_mark_reviewed_rejects_an_unknown_name():
    with pytest.raises(ValueError, match="names 0 repositories"):
        mark_reviewed(REGISTRY, "gamma", SHA_C, DAY)


def test_mark_reviewed_rejects_a_short_sha():
    with pytest.raises(ValueError, match="Not a full commit SHA"):
        mark_reviewed(REGISTRY, "beta", "abc1234", DAY)


def test_registry_rejects_an_unknown_status(tmp_path):
    path = tmp_path / "SOURCE_REPOS.toml"
    path.write_text(REGISTRY.replace('status = "watch"', 'status = "maybe"'), encoding="utf-8")
    with pytest.raises(ValueError, match="has status 'maybe'"):
        load_registry(path)


def test_owner_repo_rejects_a_url_that_is_not_a_github_repository():
    assert owner_repo("https://github.com/obra/superpowers") == "obra/superpowers"
    with pytest.raises(ValueError, match="Not a GitHub repository URL"):
        owner_repo("https://gitlab.com/obra/superpowers")
