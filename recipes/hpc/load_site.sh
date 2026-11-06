# Source from other scripts: resolve PROJECT_ROOT and the user, then load config/sites/${PROJECT_SITE}.env.
# HPC_GUIDELINES.md §1 (storage roles). PROJECT_SITE defaults to "local". USER, when unset (cron, some batch
# and container environments), is taken from `id -un`; the calling script fails here, before any side
# effect, when neither gives a name.
# shellcheck shell=bash
PROJECT_ROOT="${PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PROJECT_SITE="${PROJECT_SITE:-local}"
USER="${USER:-$(id -un 2>/dev/null || true)}"
if [[ -z "${USER}" ]]; then
    echo "ERROR: cannot determine the user name (USER is unset and id -un failed)" >&2
    exit 2
fi
site_file="${PROJECT_ROOT}/config/sites/${PROJECT_SITE}.env"
if [[ ! -f "${site_file}" ]]; then
    echo "ERROR: ${site_file} not found; copy one from recipes/hpc/sites/" >&2
    exit 2
fi
export PROJECT_ROOT PROJECT_SITE USER
# shellcheck source=/dev/null
source "${site_file}"
: "${PROJECT_WORK:?site file must define PROJECT_WORK}"
: "${PROJECT_ARCHIVE:?site file must define PROJECT_ARCHIVE}"
