#!/usr/bin/env bash
# Stage raw inputs from ARCHIVE to WORK (HPC_GUIDELINES.md §2, data flow: HPC-FLOW-001, HPC-FLOW-005).
# Usage: PROJECT_SITE=hyperion recipes/hpc/stage.sh <path under data/raw> [--dry-run]
#
# Copies ARCHIVE/data/raw/<path> to WORK/data/raw/<path> with `rsync -a --checksum`, then checks every file
# of the copy against SHA-256 checksums computed on the source before the copy; size and modification time
# are never taken as verification. <path> is relative, without '.' or '..' components, and the destination
# must resolve inside WORK/data/raw. The staged copy is a replica (CODE-PROV-001): a staged file that
# differs from the source is overwritten with the source, files present only in the staged copy are left
# in place, and the source is never modified. A verified staging appends one line to WORK/logs/staging.tsv:
# date, user, source, destination, number of files, "sha256-verified".
# --dry-run prints what would be copied and changes nothing (no directory, no log line).
set -euo pipefail
usage() { sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rel=""
dry_run=0
for arg in "$@"; do
    case "${arg}" in
        --dry-run) dry_run=1 ;;
        -h | --help) usage; exit 0 ;;
        -*) echo "ERROR: unknown option ${arg}" >&2; exit 2 ;;
        *)
            [[ -z "${rel}" ]] || { echo "ERROR: give one path under data/raw" >&2; exit 2; }
            rel="${arg}"
            ;;
    esac
done
[[ -n "${rel}" ]] || { usage >&2; exit 2; }
while [[ "${rel}" == */ ]]; do rel="${rel%/}"; done
case "/${rel}/" in
    //* | */../* | */./*)
        echo "ERROR: '${rel}' must be a relative path under data/raw, without '.' or '..'" >&2
        exit 2
        ;;
esac
# shellcheck source-path=SCRIPTDIR source=load_site.sh
source "${here}/load_site.sh"

raw_src="${PROJECT_ARCHIVE}/data/raw"
raw_dst="${PROJECT_WORK}/data/raw"
src="${raw_src}/${rel}"
dst="${raw_dst}/${rel}"
[[ -e "${src}" ]] || { echo "ERROR: ${src} does not exist" >&2; exit 2; }
src_real="$(realpath -e "${src}")"
dst_real="$(realpath -m "${dst}")"
if [[ "${src_real}" == "${dst_real}" ]]; then
    echo "source and destination are the same; nothing to stage"
    exit 0
fi
dst_root="$(realpath -m "${raw_dst}")"
if [[ "${dst_real}" != "${dst_root}/"* ]]; then
    echo "ERROR: ${dst} resolves to ${dst_real}, outside ${dst_root}" >&2
    exit 2
fi
if [[ -d "${src}" ]]; then
    from="${src}/" to="${dst}/" check_dir="${dst}"
else
    [[ ! -d "${dst}" ]] || { echo "ERROR: ${src} is a file and ${dst} a directory" >&2; exit 2; }
    from="${src}" to="${dst}" check_dir="$(dirname "${dst}")"
fi

if [[ ${dry_run} -eq 1 ]]; then
    echo "dry run: would copy ${src} to ${dst} and verify it against SHA-256 checksums of the source;"
    echo "dry run: nothing is created, copied or logged"
    if [[ -d "$(dirname "${dst}")" ]]; then
        rsync -a --checksum --itemize-changes --dry-run "${from}" "${to}"
    else
        echo "would create $(dirname "${dst}")"
        if [[ -d "${src}" ]]; then
            find "${src}" -type f -printf 'would copy %P\n'
        else
            echo "would copy ${src}"
        fi
    fi
    exit 0
fi

manifest="$(mktemp "${TMPDIR:-/tmp}/stage.XXXXXX")"
trap 'rm -f "${manifest}"' EXIT
if [[ -d "${src}" ]]; then
    (cd "${src}" && find . -type f -print0 | sort -z | xargs -0 -r sha256sum) > "${manifest}"
else
    (cd "$(dirname "${src}")" && sha256sum "./$(basename "${src}")") > "${manifest}"
fi
mkdir -p "$(dirname "${dst}")"
rsync -a --checksum --itemize-changes "${from}" "${to}"
if [[ -s "${manifest}" ]] && ! (cd "${check_dir}" && sha256sum --quiet --strict -c "${manifest}"); then
    echo "ERROR: ${dst} does not match the SHA-256 checksums of ${src}; staging not logged" >&2
    exit 6
fi
files="$(wc -l < "${manifest}")"
mkdir -p "${PROJECT_WORK}/logs"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(date -Iseconds)" "${USER}" "${src}" "${dst}" "${files}" \
    "sha256-verified" >> "${PROJECT_WORK}/logs/staging.tsv"
echo "staged ${src} -> ${dst} (${files} files verified)"
