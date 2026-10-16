"""Tests of the HPC storage scripts (HPC_GUIDELINES.md §1, storage roles; §2, data flow)."""

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


def _run(script: str, repo: Path, *args: str, site: str = "test") -> subprocess.CompletedProcess:
    env = {**os.environ, "PROJECT_SITE": site}
    env.pop("PROJECT_ROOT", None)
    return subprocess.run(
        ["bash", str(HPC / script), *args], cwd=repo, env=env, capture_output=True, text=True
    )


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
    return {"repo": repo, "work": work, "archive": archive}


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


def test_stage_copies_and_logs(project):
    repo, work, archive = project["repo"], project["work"], project["archive"]
    _run("setup_storage.sh", repo)
    (archive / "data" / "raw" / "batch1").mkdir()
    (archive / "data" / "raw" / "batch1" / "s1.fastq.gz").write_bytes(b"reads")
    result = _run("stage.sh", repo, "batch1")
    assert result.returncode == 0, result.stderr
    assert (repo / "data" / "raw" / "batch1" / "s1.fastq.gz").read_bytes() == b"reads"
    assert "batch1" in (work / "logs" / "staging.tsv").read_text()


def test_promote_writes_verified_read_only_run(project):
    repo, work, archive = project["repo"], project["work"], project["archive"]
    _run("setup_storage.sh", repo)
    run_dir = work / "results" / "20261016_test_abc1234"
    (run_dir / "tables").mkdir(parents=True)
    (run_dir / "tables" / "de.tsv").write_text("gene\tlfc\n")
    result = _run("promote.sh", repo, run_dir.name)
    assert result.returncode == 0, result.stderr
    dst = archive / "results" / run_dir.name
    manifest = (dst / "MANIFEST.sha256").read_text()
    assert "tables/de.tsv" in manifest
    assert "run: 20261016_test_abc1234" in (dst / "PROMOTION.yaml").read_text()
    assert not (dst / "tables" / "de.tsv").stat().st_mode & stat.S_IWUSR
    assert not list((archive / "results").glob("*.partial.*"))
    second = _run("promote.sh", repo, run_dir.name)
    assert second.returncode == 3
    for p in [dst, *dst.rglob("*")]:  # let pytest remove tmp_path
        p.chmod(p.stat().st_mode | stat.S_IWUSR)


def test_promote_refuses_symlinks(project):
    repo, work = project["repo"], project["work"]
    _run("setup_storage.sh", repo)
    run_dir = work / "results" / "run_with_link"
    run_dir.mkdir(parents=True)
    (run_dir / "raw").symlink_to(work / "data")
    assert _run("promote.sh", repo, run_dir.name).returncode == 5
