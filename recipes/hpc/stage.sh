#!/usr/bin/env bash
# Stage raw inputs from ARCHIVE to WORK (HPC_GUIDELINES.md §2, data flow).
# Usage: PROJECT_SITE=hyperion recipes/hpc/stage.sh <path under data/raw> [--dry-run]
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source-path=SCRIPTDIR source=load_site.sh
source "${here}/load_site.sh"
rel="${1:?usage: stage.sh <path under data/raw> [--dry-run]}"
rel="${rel%/}"
src="${PROJECT_ARCHIVE}/data/raw/${rel}"
dst="${PROJECT_WORK}/data/raw/${rel}"
[[ -e "${src}" ]] || { echo "ERROR: ${src} does not exist" >&2; exit 2; }
if [[ "$(realpath "${src}")" == "$(realpath -m "${dst}")" ]]; then
    echo "source and destination are the same; nothing to stage"
    exit 0
fi
opts=(-a --itemize-changes)
[[ "${2:-}" == "--dry-run" ]] && opts+=(--dry-run)
mkdir -p "$(dirname "${dst}")"
if [[ -d "${src}" ]]; then src="${src}/"; dst="${dst}/"; fi
rsync "${opts[@]}" "${src}" "${dst}"
mkdir -p "${PROJECT_WORK}/logs"
printf '%s\t%s\t%s\t%s\n' "$(date -Iseconds)" "${USER}" "${src}" "${dst}" >> "${PROJECT_WORK}/logs/staging.tsv"
