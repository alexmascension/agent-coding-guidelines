"""Tests of the HPC storage scripts (HPC_GUIDELINES.md §1, storage roles; §2, data flow).

The scripts run with ``USER`` removed from the environment, so the tests also check that they resolve the
user themselves (``id -un``).
"""

import hashlib
import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

HPC = Path(__file__).resolve().parents[1] / "hpc"
pytestmark = pytest.mark.skipif(
    any(shutil.which(tool) is None for tool in ("bash", "git", "rsync", "sha256sum")),
    reason="needs bash, git, rsync and sha256sum",
)


def _run(
    script: str, repo: Path, *args: str, site: str = "test", path_prefix: Path | None = None
) -> subprocess.CompletedProcess:
    env = {**os.environ, "PROJECT_SITE": site}
    for name in ("PROJECT_ROOT", "USER"):
        env.pop(name, None)
    if path_prefix is not None:
        env["PATH"] = f"{path_prefix}{os.pathsep}{env['PATH']}"
    return subprocess.run(
        ["bash", str(HPC / script), *args], cwd=repo, env=env, capture_output=True, text=True
    )


def _user() -> str:
    return subprocess.run(["id", "-un"], capture_output=True, text=True, check=True).stdout.strip()


def _fake_rsync(directory: Path, body: str) -> Path:
    """A directory holding an ``rsync`` that runs ``body`` (``$real`` is the real rsync)."""
    directory.mkdir()
    script = directory / "rsync"
    script.write_text(f'#!/usr/bin/env bash\nreal="{shutil.which("rsync")}"\n{body}\n')
    script.chmod(0o755)
    return directory


def _writable(root: Path) -> None:
    """Make a read-only promoted run writable again, so that pytest can remove tmp_path."""
    for p in [root, *root.rglob("*")]:
        p.chmod(p.stat().st_mode | stat.S_IWUSR)


@pytest.fixture
def project(tmp_path: Path) -> dict[str, Path]:
    """A git clone with a site whose WORK and ARCHIVE roots are outside it."""
    repo = tmp_path / "home" / "myrepo"
    (repo / "config" / "sites").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / ".gitignore").write_text(
        "\n".join(f"/{d}" for d in ("data", "work", "results", "logs", "images", "archive")) + "\n"
    )
    work = tmp_path / "scratch" / "user" / "myrepo"
    archive = tmp_path / "data" / "user" / "myrepo"
    (repo / "config" / "sites" / "test.env").write_text(
        f'PROJECT_WORK="{work}"\nPROJECT_ARCHIVE="{archive}"\n'
        'PROJECT_NODE="/tmp/x"\nexport PROJECT_WORK PROJECT_ARCHIVE PROJECT_NODE\n'
    )
    shutil.copy(HPC / "sites" / "local.env", repo / "config" / "sites" / "local.env")
    return {"repo": repo, "work": work, "archive": archive, "tmp": tmp_path}


@pytest.fixture
def staged(project):
    """The project after setup, with ``ARCHIVE/data/raw/batch1/s1.fastq.gz`` staged once."""
    _run("setup_storage.sh", project["repo"])
    source = project["archive"] / "data" / "raw" / "batch1" / "s1.fastq.gz"
    source.parent.mkdir()
    source.write_bytes(b"good")
    result = _run("stage.sh", project["repo"], "batch1")
    assert result.returncode == 0, result.stderr
    return {**project, "source": source, "copy": project["work"] / "data/raw/batch1/s1.fastq.gz"}


def test_setup_creates_roles_and_links(project):
    repo, work, archive = project["repo"], project["work"], project["archive"]
    result = _run("setup_storage.sh", repo)
    assert result.returncode == 0, result.stderr
    assert "not ignored" not in result.stderr
    for d in ("data/raw", "work", "results", "logs/slurm", "images", "conda", "cache"):
        assert (work / d).is_dir()
    for d in ("data/raw", "results", "images"):
        assert (archive / d).is_dir()
    for link in ("data", "work", "results", "logs", "images"):
        assert (repo / link).is_symlink()
        assert (repo / link).resolve() == (work / link).resolve()
    assert (repo / "archive").resolve() == archive.resolve()
    again = _run("setup_storage.sh", repo)
    assert again.returncode == 0 and "WARNING" not in again.stderr


def test_setup_never_replaces_existing_directory(project):
    repo = project["repo"]
    (repo / "results").mkdir()
    (repo / "results" / "keep.txt").write_text("x")
    result = _run("setup_storage.sh", repo)
    assert result.returncode == 0
    assert "exists and is not a link" in result.stderr
    assert (repo / "results" / "keep.txt").read_text() == "x"


def test_local_site_uses_real_directories(project):
    repo = project["repo"]
    result = _run("setup_storage.sh", repo, site="local")
    assert result.returncode == 0, result.stderr
    assert (repo / "data" / "raw").is_dir() and not (repo / "data").is_symlink()
    assert (repo / "archive" / "results").is_dir()


def test_stage_copies_verifies_and_logs(staged):
    repo, work = staged["repo"], staged["work"]
    assert (repo / "data" / "raw" / "batch1" / "s1.fastq.gz").read_bytes() == b"good"
    fields = (work / "logs" / "staging.tsv").read_text().rstrip("\n").split("\t")
    assert fields[1] == _user()  # USER is unset in _run: resolved with id -un
    assert fields[2].endswith("data/raw/batch1") and fields[3].endswith("data/raw/batch1")
    assert fields[4:] == ["1", "sha256-verified"]


def test_stage_stops_on_a_differing_replica_unless_overwrite(staged):
    """C09: rsync's quick check (size and mtime) would keep the wrong content; C16: replacing a replica that
    differs from the source needs --overwrite (CODE-PROV-001)."""
    source, copy = staged["source"], staged["copy"]
    copy.write_bytes(b"evil")
    times = source.stat()
    os.utime(copy, ns=(times.st_atime_ns, times.st_mtime_ns))
    log = staged["work"] / "logs" / "staging.tsv"
    refused = _run("stage.sh", staged["repo"], "batch1")
    assert refused.returncode == 7, refused.stderr
    assert "s1.fastq.gz" in refused.stderr and "--overwrite" in refused.stderr
    assert copy.read_bytes() == b"evil"
    assert len(log.read_text().splitlines()) == 1
    result = _run("stage.sh", staged["repo"], "batch1", "--overwrite")
    assert result.returncode == 0, result.stderr
    assert copy.read_bytes() == b"good"
    assert len(log.read_text().splitlines()) == 2


def test_stage_restages_an_identical_replica_without_overwrite(staged):
    result = _run("stage.sh", staged["repo"], "batch1")
    assert result.returncode == 0, result.stderr
    assert staged["copy"].read_bytes() == b"good"


def test_stage_refuses_symlinks_and_empty_sources(project):
    repo, work, archive = project["repo"], project["work"], project["archive"]
    _run("setup_storage.sh", repo)
    linked = archive / "data" / "raw" / "linked"
    linked.mkdir()
    (linked / "s1.fastq.gz").symlink_to(project["tmp"])
    assert _run("stage.sh", repo, "linked").returncode == 5
    (archive / "data" / "raw" / "empty").mkdir()
    assert _run("stage.sh", repo, "empty").returncode == 2
    assert not (work / "data" / "raw" / "linked").exists()
    assert not (work / "logs" / "staging.tsv").exists()


def test_stage_fails_when_the_copy_does_not_match_the_source(staged):
    staged["copy"].write_bytes(b"evil")
    fake = _fake_rsync(staged["tmp"] / "fake", "exit 0  # copies nothing")
    result = _run("stage.sh", staged["repo"], "batch1", "--overwrite", path_prefix=fake)
    assert result.returncode == 6
    assert "does not match the SHA-256 checksums" in result.stderr
    assert len((staged["work"] / "logs" / "staging.tsv").read_text().splitlines()) == 1


def test_stage_dry_run_creates_and_logs_nothing(project):
    repo, work, archive = project["repo"], project["work"], project["archive"]
    _run("setup_storage.sh", repo)
    for batch in ("batch1", "batch2/lane1"):
        (archive / "data" / "raw" / batch).mkdir(parents=True)
        (archive / "data" / "raw" / batch / "s1.fastq.gz").write_bytes(b"reads")
        result = _run("stage.sh", repo, batch, "--dry-run")
        assert result.returncode == 0, result.stderr
        assert "dry run" in result.stdout and "s1.fastq.gz" in result.stdout
        assert not (work / "data" / "raw" / batch).exists()
    assert not (work / "data" / "raw" / "batch2").exists()
    assert not (work / "logs" / "staging.tsv").exists()


@pytest.mark.parametrize("path", ["../outside", "/etc", "batch1/../../outside", "./batch1", ""])
def test_stage_refuses_paths_outside_data_raw(project, path):
    repo, work = project["repo"], project["work"]
    _run("setup_storage.sh", repo)
    (work / "data" / "outside").mkdir()
    result = _run("stage.sh", repo, path)
    assert result.returncode == 2
    assert not (work / "logs" / "staging.tsv").exists()


def test_stage_refuses_a_destination_that_resolves_outside_work(project):
    repo, work, archive, tmp = project["repo"], project["work"], project["archive"], project["tmp"]
    _run("setup_storage.sh", repo)
    (archive / "data" / "raw" / "batch1").mkdir()
    (archive / "data" / "raw" / "batch1" / "s1.fastq.gz").write_bytes(b"reads")
    elsewhere = tmp / "elsewhere"
    elsewhere.mkdir()
    (work / "data" / "raw" / "batch1").symlink_to(elsewhere)
    result = _run("stage.sh", repo, "batch1")
    assert result.returncode == 2
    assert list(elsewhere.iterdir()) == []


def _run_dir(work: Path, name: str = "20261016_test_abc1234") -> Path:
    run_dir = work / "results" / name
    (run_dir / "tables").mkdir(parents=True)
    (run_dir / "tables" / "de.tsv").write_text("gene\tlfc\n")
    (run_dir / "manifest.yaml").write_text("commit: abc1234\n")
    return run_dir


def test_promote_writes_verified_read_only_run(project):
    repo, work, archive = project["repo"], project["work"], project["archive"]
    _run("setup_storage.sh", repo)
    run_dir = _run_dir(work)
    result = _run("promote.sh", repo, run_dir.name)
    assert result.returncode == 0, result.stderr
    dst = archive / "results" / run_dir.name
    manifest = (dst / "MANIFEST.sha256").read_text().splitlines()
    assert [line.split("  ", 1)[1] for line in manifest] == [
        "./manifest.yaml",
        "./tables/de.tsv",
        "./PROMOTION.yaml",
    ]
    for line in manifest:  # every file of the promoted run, PROMOTION.yaml included
        digest, name = line.split("  ", 1)
        assert hashlib.sha256((dst / name).read_bytes()).hexdigest() == digest
    promotion = (dst / "PROMOTION.yaml").read_text()
    assert "run: 20261016_test_abc1234" in promotion and "files: 2" in promotion
    assert f"promoted_by: {_user()}@" in promotion  # USER is unset in _run
    assert not (dst / "tables" / "de.tsv").stat().st_mode & stat.S_IWUSR
    mode = stat.S_IMODE((dst / "PROMOTION.yaml").stat().st_mode)
    assert (
        stat.S_IMODE((dst / "MANIFEST.sha256").stat().st_mode) == mode
    )  # not 0400 from mktemp's 0600
    assert not list((archive / "results").glob("*.partial.*"))
    second = _run("promote.sh", repo, run_dir.name)
    assert second.returncode == 3
    _writable(dst)


def test_promote_checks_the_copy_against_checksums_of_the_source(project):
    """C15: a copy that differs from the source is caught, even if it is internally consistent."""
    repo, work, archive, tmp = project["repo"], project["work"], project["archive"], project["tmp"]
    _run("setup_storage.sh", repo)
    run_dir = _run_dir(work)
    fake = _fake_rsync(
        tmp / "fake",
        '"$real" "$@"\nfor last in "$@"; do :; done\nprintf x >> "${last}tables/de.tsv"',
    )
    result = _run("promote.sh", repo, run_dir.name, path_prefix=fake)
    assert result.returncode == 6, result.stderr
    assert not (archive / "results" / run_dir.name).exists()
    assert not list((archive / "results").glob("*.partial.*"))


@pytest.mark.parametrize("name", ["../outside", "a/b", ".", "..", "-x", "/abs"])
def test_promote_refuses_names_that_are_not_one_run_directory(project, name):
    """C10: nothing outside ARCHIVE/results is created or made read-only."""
    repo, work, archive = project["repo"], project["work"], project["archive"]
    _run("setup_storage.sh", repo)
    outside = work / "outside"
    outside.mkdir()
    (outside / "f.txt").write_text("x")
    result = _run("promote.sh", repo, name)
    assert result.returncode == 2
    assert not (archive / "outside").exists()
    assert (outside / "f.txt").stat().st_mode & stat.S_IWUSR
    assert list((archive / "results").iterdir()) == []


def test_promote_refuses_a_run_that_resolves_outside_results(project):
    repo, work, archive, tmp = project["repo"], project["work"], project["archive"], project["tmp"]
    _run("setup_storage.sh", repo)
    elsewhere = tmp / "elsewhere"
    elsewhere.mkdir()
    (elsewhere / "f.txt").write_text("x")
    (work / "results" / "linked").symlink_to(elsewhere)
    result = _run("promote.sh", repo, "linked")
    assert result.returncode == 2
    assert not (archive / "results" / "linked").exists()
    assert (elsewhere / "f.txt").stat().st_mode & stat.S_IWUSR


def test_promote_refuses_a_run_with_its_own_manifest(project):
    repo, work = project["repo"], project["work"]
    _run("setup_storage.sh", repo)
    run_dir = _run_dir(work)
    (run_dir / "MANIFEST.sha256").write_text("")
    assert _run("promote.sh", repo, run_dir.name).returncode == 2


def test_promote_refuses_symlinks(project):
    repo, work = project["repo"], project["work"]
    _run("setup_storage.sh", repo)
    run_dir = work / "results" / "run_with_link"
    run_dir.mkdir(parents=True)
    (run_dir / "raw").symlink_to(work / "data")
    assert _run("promote.sh", repo, run_dir.name).returncode == 5
