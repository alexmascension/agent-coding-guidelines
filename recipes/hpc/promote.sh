#!/usr/bin/env bash
# Promote one run from WORK/results to ARCHIVE/results with a manifest (HPC_GUIDELINES.md §2, data flow;
# §7, provenance on the cluster). Never overwrites an existing promoted run.
# Usage: PROJECT_SITE=hyperion recipes/hpc/promote.sh <run name>
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source-path=SCRIPTDIR source=load_site.sh
source "${here}/load_site.sh"
run="${1:?usage: promote.sh <run name>}"
run="${run%/}"
src="${PROJECT_WORK}/results/${run}"
dst="${PROJECT_ARCHIVE}/results/${run}"
tmp="${dst}.partial.$$"
[[ -d "${src}" ]] || { echo "ERROR: ${src} is not a directory" >&2; exit 2; }
[[ -e "${dst}" ]] && { echo "ERROR: ${dst} already exists; not overwriting" >&2; exit 3; }
if [[ -n "$(find "${src}" -type l -print -quit)" ]]; then
    echo "ERROR: ${src} contains symbolic links; copy their targets into the run first" >&2
    exit 5
fi
mkdir -p "$(dirname "${dst}")"
[[ -w "$(dirname "${dst}")" ]] || { echo "ERROR: $(dirname "${dst}") is not writable here" >&2; exit 4; }
trap 'rm -rf "${tmp}"' ERR

rsync -a "${src}/" "${tmp}/"
(cd "${tmp}" && find . -type f ! -name 'MANIFEST.sha256' -print0 | sort -z \
    | xargs -0 -r sha256sum > MANIFEST.sha256)
commit="$(git -C "${PROJECT_ROOT}" rev-parse HEAD 2>/dev/null || echo unknown)"
dirty="$(git -C "${PROJECT_ROOT}" status --porcelain 2>/dev/null | head -c1 | wc -c)"
{
    echo "run: ${run}"
    echo "source: ${src}"
    echo "promoted_at: $(date -Iseconds)"
    echo "promoted_by: ${USER}@$(hostname)"
    echo "repository_commit_at_promotion: ${commit}"
    echo "repository_dirty_at_promotion: $([[ ${dirty} -gt 0 ]] && echo true || echo false)"
    echo "files: $(wc -l < "${tmp}/MANIFEST.sha256")"
} > "${tmp}/PROMOTION.yaml"
(cd "${tmp}" && sha256sum --quiet -c MANIFEST.sha256)
mv "${tmp}" "${dst}"
chmod -R a-w "${dst}"
echo "promoted ${src} -> ${dst}"
