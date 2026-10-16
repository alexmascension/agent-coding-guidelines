#!/bin/bash
# Batch job template (HPC_GUIDELINES.md §4, jobs). Copy to slurm/<name>.slurm and submit from the clone
# root, where logs/ links to WORK/logs (HPC_GUIDELINES.md §1, storage roles).
# Measured on <date> (seff <jobid>): MaxRSS=<...>, Elapsed=<...>, CPU efficiency=<...>
#SBATCH --job-name=example
#SBATCH --partition=general
#SBATCH --qos=test
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=logs/slurm/%x.%j.out
#SBATCH --error=logs/slurm/%x.%j.err
#SBATCH --export=NONE
#SBATCH --no-requeue
# Arrays: add "#SBATCH --array=1-N%M" and use %x.%A_%a in --output/--error.
set -euo pipefail

# --export=NONE drops variables of the submitting shell: set the site here.
export PROJECT_SITE="${PROJECT_SITE:-hyperion}"
export PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
# shellcheck source-path=SCRIPTDIR source=load_site.sh
source "${PROJECT_ROOT}/recipes/hpc/load_site.sh"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Node-local temporary directory, removed on exit (HPC_GUIDELINES.md §6, I/O and small files).
mkdir -p "${PROJECT_NODE}"
trap 'rm -rf "${PROJECT_NODE}"' EXIT
export TMPDIR="${PROJECT_NODE}"

module load Apptainer
image="images/analysis.sif"
echo "job=${SLURM_JOB_ID} node=$(hostname) commit=$(git rev-parse HEAD) image_sha256=$(sha256sum "${image}" | cut -d' ' -f1)"

apptainer exec --bind "${PROJECT_WORK}" "${image}" \
    python -m my_project.cli run --threads "${SLURM_CPUS_PER_TASK}" --out "results/<run>/"
