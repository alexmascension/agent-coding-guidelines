# Source from other scripts: resolve PROJECT_ROOT and load config/sites/${PROJECT_SITE}.env.
# HPC_GUIDELINES.md §1 (storage roles). PROJECT_SITE defaults to "local".
# shellcheck shell=bash
PROJECT_ROOT="${PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PROJECT_SITE="${PROJECT_SITE:-local}"
site_file="${PROJECT_ROOT}/config/sites/${PROJECT_SITE}.env"
if [[ ! -f "${site_file}" ]]; then
    echo "ERROR: ${site_file} not found; copy one from recipes/hpc/sites/" >&2
    exit 2
fi
export PROJECT_ROOT PROJECT_SITE
# shellcheck source=/dev/null
source "${site_file}"
: "${PROJECT_WORK:?site file must define PROJECT_WORK}"
: "${PROJECT_ARCHIVE:?site file must define PROJECT_ARCHIVE}"
