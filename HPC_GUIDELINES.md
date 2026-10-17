# HPC guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 1.1.0 (2026-10-17)
>
> Key words **MUST**, **SHOULD**, **MAY** follow RFC 2119. Each rule is followed by its rationale: ***Why***
> (mechanism), ***Trade-off*** (cost) and ***Deviate when*** (justified exceptions), where they apply.
> A rule without a key word is a default with the force of SHOULD. Rules carry stable identifiers (e.g.
> `HPC-STOR-001`); overrides refer to them (`CODING_GUIDELINES.md` §0, rule identifiers).

## 0. How to use this document

* Audience: coding agents and people who run analyses on shared clusters with a batch scheduler (Slurm).
* Copy with `CODING_GUIDELINES.md`. This file is site-independent: it defines **roles** of storage and
  conventions for jobs. The facts of one cluster (mount points, quotas, partitions, QoS, modules) are in a **site
  profile**, `hpc_sites/<SITE>.md`, copied into the project (`guidelines/hpc_sites/`, `CODE-USE-001`) and referenced from the *Stack* section of
  `PROJECT_CONTEXT.md`. Site facts change: profiles carry the date they were checked against the centre's
  documentation, and the centre's documentation wins when they disagree.
* Scope: batch analysis (pipelines, R/Python jobs, GPU jobs) of research projects. Not covered: MPI development,
  cluster administration.

## 1. Storage roles and layout

* `HPC-STOR-001` **Every project has the same relative path on every storage area**: `<root>/$USER/<repo>`, where `<repo>` is the
  repository name and `<root>` is the mount point of each role. Roles:

  | Role | Holds | Typical properties |
  |---|---|---|
  | `CODE` | The git clone: code, configuration, small documents | Home; small quota, many-file limit |
  | `WORK` | Everything jobs read and write: staged inputs, intermediates, results, logs, environments, caches | Parallel scratch filesystem; large quota; not permanent, often not backed up |
  | `ARCHIVE` | Raw inputs and promoted final results kept long term | Long-term storage; may be read-only from compute nodes |
  | `NODE` | Temporary files of one job: `<root>/$USER/$SLURM_JOB_ID` | Node-local disk; fastest; deleted (or orphaned) when the job ends |

  * *Why:* with one naming rule, every path of a project is predictable from the role and the repository name, for
    people and for agents; scripts need only the role roots, and two projects never collide.
  * *Deviate when:* a role does not exist on a site (e.g. no `ARCHIVE`); the site profile states the substitute.
* `HPC-STOR-002` **Fixed subdirectories inside each role root** (`<root>/$USER/<repo>/...`): `WORK`: `data/` (staged inputs, same
  relative paths as in `ARCHIVE`, e.g. `data/raw/...`), `work/` (pipeline intermediates), `results/`, `logs/`,
  `images/` (container images), `conda/` (conda environments, if any), `cache/` (package and image
  caches); `ARCHIVE`: `data/raw/` (original inputs, read-only), `results/` (promoted runs), `images/` (images of
  promoted results).
* `HPC-STOR-003` **Code finds the roles through variables, never through absolute paths** (`CODE-CONF-002`): `PROJECT_WORK`,
  `PROJECT_ARCHIVE`, `PROJECT_NODE`, and `PROJECT_ROOT` for the clone (`CODE-CONF-003`). A tracked file per site,
  `config/sites/<site>.env`, defines them from the role roots and `$USER`; `PROJECT_SITE` selects the file, and the
  default site `local` places every role inside the repository (gitignored directories), so a fresh clone runs on a
  laptop (`CODE-CONF-001`).
* `HPC-STOR-004` **The repository exposes the roles as gitignored symbolic links** created by a tracked script
  (`recipes/hpc/setup_storage.sh`): `data/`, `work/`, `results/`, `logs/` and `images/` in the clone point to the
  same directories of `WORK`, and `archive/` points to the `ARCHIVE` root. Code uses project-relative paths
  (`data/raw/...`, `results/<run>/...`) and works unchanged on the cluster and locally; only staging and promotion
  touch `archive/` (see §2, data flow).
  * *Why:* relative paths keep notebooks and scripts identical everywhere; the links make the split across
    filesystems visible with `ls -l`.
* `HPC-STOR-005` **Code, environments and caches never go to `CODE` beyond source and configuration.** Conda environments, pip
  and container caches, Nextflow work directories and job output go to `WORK`.
  * *Why:* home areas have small quotas and file-count limits; one conda environment holds tens of thousands of
    files.

## 2. Data flow, staging and backups

* `HPC-FLOW-001` **Raw inputs live once in `ARCHIVE/data/raw` and are read-only** (`CODE-PROV-001`). Jobs SHOULD read staged copies
  in `WORK/data/` when `ARCHIVE` is not designed for job I/O (e.g. an NFS filesystem, or read-only from compute
  nodes); staging is an explicit, logged `rsync` step.
  * *Trade-off:* staged copies count against the `WORK` quota; stage what the next runs need and remove it after.
  * *Deviate when:* large read-only references (e.g. aligner indexes) MAY be read directly from `ARCHIVE` when load
    times measured under the intended concurrency are acceptable; record the measurement.
* `HPC-FLOW-002` **Jobs write only to `WORK` (and `NODE`).** Final results reach `ARCHIVE` through an explicit **promotion** step
  (`recipes/hpc/promote.sh`), run where `ARCHIVE` is writable (e.g. a login node), that copies one run directory
  and writes a manifest (file list with checksums, commit, job IDs, date).
  * *Why:* results in `ARCHIVE` are then complete runs with provenance, never half-written outputs of a job that
    failed; and the rule works on sites where `ARCHIVE` is read-only from compute nodes.
* `HPC-FLOW-003` **Assume nothing on the cluster is backed up unless the site profile says so.** Code is safe only when pushed to
  the git remote; raw data and promoted results that cannot be regenerated MUST have a copy outside the cluster
  (institutional storage, repository deposit), recorded in `PROJECT_CONTEXT.md`.
* `HPC-FLOW-004` **`WORK` is not storage.** Remove intermediates of finished runs (or archive them as one file); scratch areas
  degrade for everyone when full and may be purged by the centre.
* `HPC-FLOW-005` **Transfers use `rsync`** (`-a`, checksums for verification, `--dry-run` before any `--delete`); large or many
  files are archived (`tar`) before transfer. Downloads from public archives follow `CODING_GUIDELINES.md` §13
  (downloads from public archives).

### Datasets larger than the work area

* `HPC-BATCH-001` **Process in batches sized to a space budget.** Each batch fetches or stages its inputs, runs the pipeline,
  harvests and validates its results, and deletes its inputs and intermediates before the next batch takes the space.
  The budget is the project's share of the `WORK` quota (shared with the user's other projects) minus a free margin;
  a batch waits, for a bounded time, when the budget is exhausted.
* `HPC-BATCH-002` **Calibrate the footprint on a pilot**: peak `WORK` usage divided by input size, measured on a few batches;
  record it in `PROJECT_CONTEXT.md` and size batches from it.
* `HPC-BATCH-003` **Stage markers make a batch restartable**: one marker file per completed stage, written after the stage's outputs
  are complete; a relaunch skips completed stages, and workflow managers resume (`-resume`) with their work
  directory inside the batch.
* `HPC-BATCH-004` **MUST validate before deleting**: count and check the expected outputs of every sample of the batch (and
  promote or harvest them) before any input or intermediate is deleted; on any mismatch, delete nothing and fail.
* `HPC-BATCH-005` **Extract everything that needs an intermediate before deleting it** (e.g. genotype fingerprints and region BAMs
  before BAMs are deleted; `SAMPLE_IDENTITY_GUIDELINES.md` §4, genotype fingerprints).
* `HPC-BATCH-006` Concurrent batches MAY be capped without a controller process with `--dependency=singleton`: one job name per
  slot, so at most one batch per name runs at a time.

## 3. Software environments

* `HPC-ENV-001` **Default: one Apptainer (Singularity) image per environment, built from tracked definition files**
  (`envs/<name>.def` and its `environment.yml` or lock file, tracked in the repository), with the built `.sif` in
  `WORK/images/` and, for promoted results, in `ARCHIVE/images/`.
  * *Why:* a conda environment is tens of thousands of small files, slow to activate on parallel filesystems and
    heavy on file-count quotas; an image is one file, and the definition file reproduces it.
  * *Trade-off:* rebuilding the image for every dependency change; building may need a node or machine where
    builds are allowed.
  * *Deviate when:* small or short-lived environments MAY be conda environments created with `--prefix` in
    `WORK/conda/` from a tracked `environment.yml`; the site's own modules MAY be used
    for tools the centre provides (record module names and versions).
* `HPC-ENV-002` **Record the image digest (`sha256sum` of the `.sif`) with every result** (see §7, provenance on the
  cluster); a tag or file name is not enough.
* `HPC-ENV-003` **Redirect every cache to `WORK/cache/`** before building or pulling: `APPTAINER_CACHEDIR`, conda `pkgs_dirs`,
  `PIP_CACHE_DIR`, Nextflow `NXF_HOME` and `NXF_APPTAINER_CACHEDIR` (or `apptainer.cacheDir`). Caches written by
  jobs MUST be writable from compute nodes.
* `HPC-ENV-004` GPU jobs run images with `--nv`; the CUDA version of the image MUST be supported by the node driver (check on a
  test job).

## 4. Jobs

* `HPC-JOB-001` **Batch scripts are tracked in the repository** (`slurm/` or the workflow manager configuration), take paths and
  parameters as arguments, and start from the template `recipes/hpc/job.slurm`:
  * `#SBATCH --export=NONE` and explicit loading of modules or images inside the script;
  * `set -Eeuo pipefail` and an `ERR` trap that reports the failing line (`CODE-LANG-004`);
  * `--output`/`--error` to `WORK/logs/` with job name and ID (`%x.%j`, arrays `%x.%A_%a`);
  * `--no-requeue` unless the job resumes from checkpoints;
  * thread counts from the allocation (`OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK`, and the tool's own thread option).
  * *Why:* propagating the login environment makes jobs depend on whatever was loaded at submission; automatic
    requeueing of a job that does not resume restarts it from the beginning and overwrites its outputs.
* `HPC-JOB-002` **Resources are requested from measurements**: run once with generous limits on a small input, read `seff` and
  `sacct` (`MaxRSS`, `Elapsed`, CPU efficiency) of completed jobs, then request about 1.2–1.5 times the measured
  peak memory and time. Record the measured values in the job script comments.
  * *Why:* oversized requests wait longer in the queue and waste fair share; undersized ones fail at the end.
* `HPC-JOB-003` **Test before production**: a short run on the test QoS or a reduced input, then the full job.
* `HPC-JOB-004` **Independent tasks use job arrays** with a throttle (`--array=1-500%50`) and one task per sample, not one job per
  sample submitted in a loop.
  * *Why:* arrays respect submission limits and are cancelled, monitored and accounted as one job.
* `HPC-JOB-005` **Jobs are idempotent and restartable**: each task writes to a temporary name and renames on success, skips work
  whose output already exists and is complete, and never overwrites existing outputs without `--overwrite`
  (`CODE-PROV-005`).
* `HPC-JOB-006` **Preemptible partitions only for jobs that tolerate being killed at any time**: checkpointing jobs or short
  idempotent array tasks.
* `HPC-JOB-007` **Interactive work runs in an allocation** (`salloc`/`srun --pty bash`), never on login nodes.

## 5. Workflow managers

* `HPC-WF-001` **Nextflow, Snakemake or similar submit each step as a job**; the controller process is light but long-lived:
  run it as its own long, low-resource job, or in `tmux`/`screen` on a login node if the site allows it.
* `HPC-WF-002` **Work directory, caches and outputs in `WORK`** (`-w $PROJECT_WORK/work/<pipeline>`); outputs then promoted as in
  §2 (data flow). Use `-resume` after failures instead of restarting from zero.
* `HPC-WF-003` **Pin pipeline and engine versions** (`-r <release>` for nf-core pipelines, the module version of Nextflow) and keep the
  site configuration (executor, partition, QoS, container settings) in a tracked config file, separate from the
  pipeline parameters.
* `HPC-WF-005` **Pipeline parameters go in a tracked parameters file** (`-params-file params.yaml`), not on the command line.
  * *Why:* the file is the record of the run; and with Nextflow 26, boolean and numeric parameters given on the
    command line reached nf-core schema validation as strings in at least one project, while a YAML file keeps
    their types.
* `HPC-WF-006` **With Apptainer, check that every process uses a prebuilt image.** Some nf-core modules recognise only the
  `singularity` engine and, under `-profile apptainer`, fall back to `docker://` images that are built from layers on
  the fly, which can fail on parallel filesystems; `-profile singularity` (Apptainer provides the `singularity`
  command) selects the prebuilt `.sif` images. Pin individual processes to an image in the site configuration when
  needed.
* `HPC-WF-004` Clean the work directory of completed runs (`nextflow clean`) once results are promoted.

## 6. I/O and small files

* `HPC-IO-001` **Avoid many small files on parallel filesystems**: keep fewer than a few thousand files per directory; pack
  per-sample outputs into archives or container formats (HDF5, zarr, parquet); use `NODE` for workloads that create
  many temporary files and copy back only the results.
* `HPC-IO-002` **`NODE` directories are created and removed by the job** (`mkdir -p "$PROJECT_NODE"` at the start, cleanup in a
  `trap ... EXIT`), with the job ID in the name; a killed job may leave them behind, so the name must identify them.
* `HPC-IO-003` Read large inputs sequentially and once per job; avoid many jobs reading the same file from a non-parallel
  filesystem at the same time (stage it first; §2, data flow).

## 7. Provenance on the cluster

* `HPC-PROV-001` **Every run directory has a manifest** written by the job or the workflow manager: commit of the repository (and
  whether the tree was clean), parameters, image digest or module versions, Slurm job IDs (and array indices),
  node names, start and end time, and the requested and used resources (`sacct`). The promotion step of
  §2 (data flow) copies it with the results.
  * *Why:* a result on a shared cluster depends on the code, the environment and the node it ran on; without the
    manifest none of them can be reconstructed after the scratch area is cleaned.
* `HPC-PROV-002` Runs are named `<YYYYMMDD>_<short description>_<short commit>`; the same name is used in `WORK/results/` and in
  `ARCHIVE/results/`.

## 8. Agents on shared clusters

* `HPC-AGENT-001` **MUST NOT run computations on login nodes.** Inspect files, edit, submit, and run light commands only; anything
  that uses noticeable CPU, memory or I/O for more than a few seconds goes into a job. Centres terminate such
  processes and may suspend the account.
* `HPC-AGENT-002` **MUST ask before** submitting more than a handful of jobs or any job array larger than the project's agreed
  threshold, before using a partition or QoS not listed in `PROJECT_CONTEXT.md`, and before deleting or
  overwriting anything in `WORK` or `ARCHIVE` (`CODE-AGENT-004`).
* `HPC-AGENT-003` **Cancel only own jobs, by job ID**; never `scancel -u $USER` without confirmation.
* `HPC-AGENT-004` Check quotas (size and file count) before jobs that create large outputs, and report usage when it exceeds about
  80% of a quota.
* `HPC-AGENT-005` Report job IDs, the expected run time and how to monitor them (`squeue`, `sacct`) after submitting.

## 9. What to record in `PROJECT_CONTEXT.md`

1. `HPC-REP-001` Site profile used and its date; the values of the role roots; partitions and QoS the project may use, including
   any granted to the group.
2. `HPC-REP-002` Image names, definition files and digests; modules used.
3. `HPC-REP-003` Staging and promotion rules of the project (which inputs are staged, which outputs are promoted).
4. `HPC-REP-004` Off-cluster copies of raw data and promoted results.
5. `HPC-REP-005` Thresholds for agents (maximum jobs or array size without confirmation).

## 10. Implementation

| Task | File |
|---|---|
| Role roots of a site (copy to `config/sites/<site>.env`) | `recipes/hpc/sites/local.env`, `recipes/hpc/sites/hyperion.env` |
| Resolve the project root and load the site file (sourced by the scripts below) | `recipes/hpc/load_site.sh` |
| Create role directories and repository links (§1) | `recipes/hpc/setup_storage.sh` |
| Stage raw inputs from `ARCHIVE` to `WORK`, with a log (§2) | `recipes/hpc/stage.sh` |
| Promote a run from `WORK` to `ARCHIVE` with checksums and a manifest (§2, §7) | `recipes/hpc/promote.sh` |
| Batch job template with `ERR` trap and node-local directory (§4, §6) | `recipes/hpc/job.slurm` |
| Site facts | `hpc_sites/<SITE>.md` (e.g. `hpc_sites/DIPC_HYPERION.md`) |

Tests: `recipes/tests/test_hpc_scripts.py` (temporary directories as role roots; no scheduler needed). The job
template is checked with `shellcheck` only.
