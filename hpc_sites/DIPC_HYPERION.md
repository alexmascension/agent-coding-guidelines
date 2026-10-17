# Site profile: DIPC Hyperion

> Companion of `HPC_GUIDELINES.md`. Checked against the DIPC Supercomputing Center documentation
> (<https://scc.dipc.org/docs/systems/hyperion/>, <https://scc.dipc.org/docs/jobs/slurm/>,
> <https://scc.dipc.org/docs/jobs/gpus/>, <https://scc.dipc.org/docs/getting-started/storage/>) on **2026-10-16**; §8 adds observations
> made by a project on the cluster (2026-10-03 to 2026-10-17). Facts, not rules: when this file and the centre's documentation disagree, the documentation wins; update the file
> and its date. Support: `support-hpc@dipc.org`.

## 1. Access and login nodes

* `ssh <user>@hyperion.sw.ehu.es` (login nodes `hyperion-01`, `hyperion-02`).
* Login nodes are for file management, editing and compilation. **Resource-intensive processes on login nodes are
  terminated and the account may be suspended** (`HPC-AGENT-001`). Interactive work:
  `salloc --qos=test ...` or `srun --qos=test --pty bash`.
* Long-running controllers (Nextflow) run inside `tmux` or `screen`.

## 2. Storage roles

| Role (`HPC_GUIDELINES.md` §1) | Hyperion path | Quota (size / files) | Notes |
|---|---|---|---|
| `CODE` | `/dipc/$USER/<repo>` (home) | 50 GB / 300K | Shared by login and compute nodes. Configuration, dotfiles, small data; no environments or job output. |
| `WORK` | `/scratch/$USER/<repo>` | 2 TB / 10M | BeeGFS over InfiniBand; all job I/O; submit jobs from here. Not permanent; performance degrades for everyone above 80% occupancy. `/scratch/$USER` is created with the account. |
| `ARCHIVE` | `/data/$USER/<repo>` | 10 TB / 20M | NFS over the 100 GbE management network; **read-only on compute nodes, writable only from login nodes**. For data not actively used by running jobs. One-day grace period up to 10.5 TB / 20.5M files, with an email alert. |
| `NODE` | `/lscratch/$USER/$SLURM_JOB_ID` | node disk, 750 GB – 14 TB depending on the node | Node-local; not created automatically; create at job start and remove at the end. Absent on login nodes; check that it exists (§8). |

* **No filesystem on Hyperion is backed up.** Code on the git remote and off-cluster copies of raw data and
  promoted results are the only protection (`HPC_GUIDELINES.md` §2).
* Consequences for the data flow of `HPC_GUIDELINES.md` §2: jobs cannot write to `/data`; **staging** (`/data` →
  `/scratch`) and **promotion** (`/scratch` → `/data`) run on a login node with `rsync`; jobs read staged copies in
  `/scratch`, not `/data`.
* Quotas (size and file count): `getquota` (all), `getquota-home`, `getquota-scratch`, `getquota-data`. A file count at
  the limit blocks new files even with free space. Temporary quota increases are granted on request.
* Avoid more than a few thousand files per directory on `/scratch`; pack them (`tar`) or use `/lscratch`.
* Sharing with other users: POSIX ACLs (`setfacl`) on `/scratch` or `/data` directories.

## 3. Partitions and QoS

These are the **general** partitions and QoS available to every account. Groups can be granted further partitions,
QoS or temporary limits (on request to support); record them in `PROJECT_CONTEXT.md`. Check the live values with
`sacctmgr show qos` and `sinfo`.

| Partition | Nodes | Notes |
|---|---|---|
| `general` (default) | All public nodes | |
| `preemption` | Group-owned nodes open to everyone | Jobs can be killed and requeued when the owners need the nodes; use a QoS of the general partition. Some nodes have CPU architectures (Zen 3, Zen 4) outside the general set (Cascade Lake, Ice Lake, Emerald Rapids); centrally installed software is not guaranteed there. |
| `preemption-gpu` | A100 SXM4 nodes | Small GPU jobs (< 4 GPUs), preemptible; use its default QoS only. |

QoS for `general` and `preemption`:

| QoS | Priority | Max walltime | Max nodes/user | Max running jobs/user | Max submitted | Max TRES |
|---|---|---|---|---|---|---|
| `regular` (default) | 200 | 1 day | 60 | 180 | | |
| `test` | 1000 | 10 min | 2 | 2 | 2 | |
| `long` | 200 | 2 days | 25 | 40 | | |
| `xlong` | 200 | 8 days | 20 | 20 | 200 | |
| `serial` | 200 | 2 days | | 1000 | 2000 | CPUs=1, GPUs=0, nodes=1 |

QoS for `preemption-gpu`: `preemption-gpu` (default), priority 200, 8 days, GPUs=3, nodes=1.

* **Global limit: 1000 submitted (queued + running) jobs across all QoS**; `serial` allows up to 2000.
* Job arrays: at most 10,000 tasks per array and indices below 4,000,000; throttle with `%N`.
* Priority combines fair share (recent usage relative to the assigned share), age in queue and QoS priority;
  backfill can start short jobs earlier, so accurate `--time` requests start sooner.
* Typical choices: tests and resource measurements → `test`; one sample per task, single-threaded → `serial`;
  multi-threaded per-sample tasks (alignment, quantification) → `regular`; jobs longer than 1 day → `long`/`xlong`
  only after checking they cannot be split.

## 4. Job script conventions on Hyperion

* `#SBATCH --export=NONE` and modules loaded inside the script; `#SBATCH --no-requeue` unless the job resumes from
  checkpoints (after a node failure Slurm requeues by default and the restarted job overwrites outputs).
* Submit from and write output to `/scratch/$USER/<repo>`.
* Measure resources with `seff <jobid>`, `sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,ReqMem,AllocCPUS,State`,
  and `sstat` for running jobs.
* Node availability: `cluster-resources` (filters: `cluster-resources --help`).

## 5. GPUs

* Request with `--gres=gpu:N` and select the type with `--constraint`:

  | GPU | Memory | GPUs/node | Constraint | Partition |
  |---|---|---|---|---|
  | A100 SXM4 (NVLink) | 80 GB | 8 | `a100-sxm4` | `general` (≥ 4 GPUs per node), `preemption-gpu` (< 4) |
  | A100 PCIe | 80 GB | 8 | `a100-pcie` | `general` |
  | RTX 3090 | 24 GB | 2 | `rtx3090` | `general` (also `preemption`) |
  | RTX A6000 | 48 GB | 1 | `a6000` | `general` |
  | A30 | 24 GB | 2 | `a30` | `preemption` |
  | RTX PRO 6000 | 96 GB | 1 | `pro6000` | `preemption` |
  | RTX 4090 | 24 GB | 1 | `rtx4090` | `preemption` |

* Jobs requesting fewer than 4 GPUs cannot use A100 SXM4 nodes on `general`; use `a100-pcie` or `preemption-gpu`.
* CUDA from modules (e.g. `module load CUDA/12.4.0`) or inside an image run with `apptainer exec --nv`.

## 6. Software

* Modules: `module avail`, `module spider <name>`, `module load <name>/<version>`; record exact versions.
* Containers: `module load Apptainer`; set `APPTAINER_CACHEDIR=/scratch/$USER/.apptainer-cache` (or the project's
  `WORK/cache/apptainer`); the centre recommends packaging conda environments inside images (`HPC-ENV-001`).
* Conda: `module load Miniforge3`; environments with `--prefix` under `/scratch/$USER/<repo>/conda/` (an environment under
  `/data` must be created from a login node, because `/data` is read-only on compute nodes).
* Nextflow (see also §8): `module load Nextflow/25.04.6`; `executor = 'slurm'`; `queue` sets the partition and
  `clusterOptions = '--qos=...'` the QoS; work directory `-w /scratch/$USER/<repo>/work/nf`; Apptainer loaded in the
  session and in each task, `apptainer.cacheDir` on `/scratch`; clean finished runs with `nextflow clean -f`.

## 7. Site file

`recipes/hpc/sites/hyperion.env` (copy to `config/sites/hyperion.env`) defines the role roots above and redirects
the Apptainer, pip and conda caches to `/scratch/$USER/<repo>/cache/`. Run `recipes/hpc/setup_storage.sh` with
`PROJECT_SITE=hyperion` **on a login node**: it creates `/data/$USER/<repo>`, which compute nodes cannot write.
Staging (`recipes/hpc/stage.sh`) and promotion (`recipes/hpc/promote.sh`) also run on a login node.
Submitting from the clone in `/dipc` is compatible with the centre's advice to submit from `/scratch` as long as
every output path of the job resolves to `/scratch` through the links of `HPC-STOR-004`.

## 8. Observed on the cluster

Observations of one project (bulk RNA-seq with nf-core/rnaseq 3.27.0), not statements of the centre. Re-check them
before relying on them, and remove them when the documentation covers them.

* **Nextflow version:** the module `Nextflow/25.04.6` was too old for nf-core/rnaseq 3.27.0 (`Plugin nf-schema with
  version @2.7.2 does not exist`). A launcher pinned with `NXF_VER` (26.04.7, then 26.10.0) in `WORK` worked, after
  `module load Java` (Java 17 or later; 21 available).
* **`NXF_HOME`** under `/data` made Nextflow fail at start-up on compute nodes (read-only); keep it in
  `WORK/cache/nextflow` (`recipes/hpc/sites/hyperion.env`).
* **nf-core with Apptainer:** use `-profile singularity` (after `module load Apptainer`, `singularity` is provided by
  Apptainer). With `-profile apptainer`, some modules (e.g. `DESEQ2_QC` of nf-core/rnaseq 3.27.0) fell back to
  `docker://` images, and building them from layers on BeeGFS failed with `unpriv.link: operation not permitted`.
* **Reading `/data` from jobs:** a 29 GB STAR index loaded in a median of 12 s from `/data` versus 13 s from `/scratch`
  (20 tasks, alternating; one outlier of 62 s from `/data`). Reading large read-only references directly from
  `/data` was viable for one batch at a time; concurrency was not measured (`HPC-FLOW-001`, *Deviate when*).
* **Quota from scripts:** `beegfs-ctl --getquota --uid … --csv` returned nothing for a regular user; use
  `getquota-scratch` or `getquota`.
* **Network:** compute nodes reach the internet through a different address than login nodes, and an archive
  (ENA) served broken responses for some files to compute nodes only, for at least 15–20 minutes. Test downloads
  on a compute node (`CODING_GUIDELINES.md` §13, downloads from public archives).
* **`/lscratch`:** absent on login nodes and not found on some compute nodes, although the documentation lists a
  size for most nodes; on `hyperion-294` (14 TB NVMe) `/tmp` points to `/lscratch/tmp`. Scripts check that it exists
  and fall back to `WORK` (`recipes/hpc/sites/hyperion.env`).
* **Group partitions** can require the account and the QoS together (`sbatch -p <partition> -A <account>
  --qos=<qos>`); `scontrol show partition <partition>` lists `AllowAccounts` and `AllowQos`. With a QoS that is not
  allowed, `srun` printed `The requested QoS is not allowed in the used partition` and changed it. Record such
  partitions in `PROJECT_CONTEXT.md` (`HPC-REP-001`).
* **Home path:** `/dipc/$USER` resolves to `/home/$USER` (`readlink -f`), and Slurm records `WorkDir=/home/…` for jobs
  submitted from the clone; it is the same clone, not a second one.
