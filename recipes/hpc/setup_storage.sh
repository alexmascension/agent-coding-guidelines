#!/usr/bin/env bash
# Create the role directories of a project and link them into the clone.
# HPC_GUIDELINES.md §1 (storage roles). Run from the clone; on sites where ARCHIVE is read-only on
# compute nodes, run on a login node. Idempotent; never replaces an existing file or directory.
# Usage: PROJECT_SITE=hyperion recipes/hpc/setup_storage.sh [--dry-run]
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source-path=SCRIPTDIR source=load_site.sh
source "${here}/load_site.sh"
dry_run=0
[[ "${1:-}" == "--dry-run" ]] && dry_run=1

run() {
    echo "+ $*"
    if [[ ${dry_run} -eq 0 ]]; then "$@"; fi
}

work_dirs=(data/raw work results logs/slurm images conda cache)
archive_dirs=(data/raw results images)
links=(data work results logs images)

for d in "${work_dirs[@]}"; do run mkdir -p "${PROJECT_WORK}/${d}"; done
for d in "${archive_dirs[@]}"; do run mkdir -p "${PROJECT_ARCHIVE}/${d}"; done

link() {  # link <path in clone> <target>
    local path="${PROJECT_ROOT}/$1" target="$2"
    if [[ "$(realpath -m "${path}")" == "$(realpath -m "${target}")" ]]; then
        return 0  # role root is the clone itself (local site): nothing to link
    fi
    if [[ -L "${path}" ]]; then
        if [[ "$(readlink "${path}")" != "${target}" ]]; then
            echo "WARNING: ${path} links to $(readlink "${path}"), expected ${target}; left unchanged" >&2
        fi
    elif [[ -e "${path}" ]]; then
        echo "WARNING: ${path} exists and is not a link; left unchanged" >&2
    else
        run ln -s "${target}" "${path}"
    fi
}

for l in "${links[@]}"; do link "${l}" "${PROJECT_WORK}/${l}"; done
link archive "${PROJECT_ARCHIVE}"

for l in "${links[@]}" archive; do
    if ! git -C "${PROJECT_ROOT}" check-ignore -q "${l}"; then
        echo "WARNING: '${l}' is not ignored by git; add '/${l}' to .gitignore" >&2
    fi
done
echo "site=${PROJECT_SITE} work=${PROJECT_WORK} archive=${PROJECT_ARCHIVE}"
