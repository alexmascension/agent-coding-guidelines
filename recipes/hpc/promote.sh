#!/usr/bin/env bash
# Promote one run from WORK/results to ARCHIVE/results with a manifest (HPC_GUIDELINES.md §2, data flow;
# §7, provenance on the cluster). Never overwrites an existing promoted run.
# Usage: PROJECT_SITE=hyperion recipes/hpc/promote.sh <run name>
#
# <run name> is one directory of WORK/results (letters, digits, '.', '_', '-'; no '/'), and source and
# destination must resolve directly under WORK/results and ARCHIVE/results. SHA-256 checksums are computed
# on the run in WORK before the copy (MANIFEST.sha256) and the copy is checked against them, so a file that
# was copied wrongly, or changed in WORK during the copy, stops the promotion. PROMOTION.yaml is written
# next and its checksum appended to MANIFEST.sha256: `sha256sum -c MANIFEST.sha256` in the promoted
# directory checks every file. The copy goes to a temporary directory next to the destination, is renamed
# once verified, and is made read-only.
set -euo pipefail
usage() { sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
case "${1:-}" in -h | --help) usage; exit 0 ;; esac
[[ $# -eq 1 && -n "${1}" ]] || { usage >&2; exit 2; }
run="${1%/}"
if [[ ! "${run}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "ERROR: '${run}' is not one directory name (letters, digits, '.', '_', '-'; no '/')" >&2
    exit 2
fi
# shellcheck source-path=SCRIPTDIR source=load_site.sh
source "${here}/load_site.sh"

src="${PROJECT_WORK}/results/${run}"
dst_parent="${PROJECT_ARCHIVE}/results"
dst="${dst_parent}/${run}"
[[ -d "${src}" ]] || { echo "ERROR: ${src} is not a directory" >&2; exit 2; }
src_root="$(realpath -e "${PROJECT_WORK}/results")"
src_real="$(realpath -e "${src}")"
if [[ "$(dirname "${src_real}")" != "${src_root}" ]]; then
    echo "ERROR: ${src} resolves to ${src_real}, not a directory of ${src_root}" >&2
    exit 2
fi
dst_real="$(realpath -m "${dst}")"
if [[ "$(dirname "${dst_real}")" != "$(realpath -m "${dst_parent}")" ]]; then
    echo "ERROR: ${dst} resolves to ${dst_real}, outside ${dst_parent}" >&2
    exit 2
fi
if [[ -e "${dst}" || -L "${dst}" ]]; then
    echo "ERROR: ${dst} already exists; not overwriting" >&2
    exit 3
fi
if [[ -n "$(find "${src}" -type l -print -quit)" ]]; then
    echo "ERROR: ${src} contains symbolic links; copy their targets into the run first" >&2
    exit 5
fi
for reserved in MANIFEST.sha256 PROMOTION.yaml; do
    if [[ -e "${src}/${reserved}" ]]; then
        echo "ERROR: ${src}/${reserved} exists; promotion writes that file" >&2
        exit 2
    fi
done
mkdir -p "${dst_parent}"
[[ -w "${dst_parent}" ]] || { echo "ERROR: ${dst_parent} is not writable here" >&2; exit 4; }

tmp="${dst}.partial.$$"
manifest="$(mktemp "${TMPDIR:-/tmp}/promote.XXXXXX")"
trap 'rm -rf "${tmp}" "${manifest}"' EXIT
(cd "${src}" && find . -type f -print0 | sort -z | xargs -0 -r sha256sum) > "${manifest}"
files="$(wc -l < "${manifest}")"
[[ "${files}" -gt 0 ]] || { echo "ERROR: ${src} contains no files" >&2; exit 2; }

rsync -a "${src}/" "${tmp}/"
cp "${manifest}" "${tmp}/MANIFEST.sha256"
copied="$(cd "${tmp}" && find . -type f ! -path ./MANIFEST.sha256 -printf . | wc -c)"
if [[ "${copied}" -ne "${files}" ]] \
    || ! (cd "${tmp}" && sha256sum --quiet --strict -c MANIFEST.sha256); then
    echo "ERROR: the copy does not match the checksums computed on ${src} (${copied} files copied," \
        "${files} in the manifest); nothing promoted" >&2
    exit 6
fi
commit="$(git -C "${PROJECT_ROOT}" rev-parse HEAD 2>/dev/null || echo unknown)"
dirty="$(git -C "${PROJECT_ROOT}" status --porcelain 2>/dev/null | head -c1 | wc -c)"
{
    echo "run: ${run}"
    echo "source: ${src}"
    echo "promoted_at: $(date -Iseconds)"
    echo "promoted_by: ${USER}@$(hostname)"
    echo "repository_commit_at_promotion: ${commit}"
    echo "repository_dirty_at_promotion: $([[ ${dirty} -gt 0 ]] && echo true || echo false)"
    echo "files: ${files}"
    echo "manifest: MANIFEST.sha256, computed on the source before the copy; PROMOTION.yaml appended"
} > "${tmp}/PROMOTION.yaml"
(cd "${tmp}" && sha256sum ./PROMOTION.yaml >> MANIFEST.sha256 && sha256sum --quiet --strict -c MANIFEST.sha256)
mv -T "${tmp}" "${dst}"
chmod -R a-w "${dst}"
echo "promoted ${src} -> ${dst} (${files} files verified)"
