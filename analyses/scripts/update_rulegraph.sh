#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
graph_tmp=$(mktemp analyses/docs/.rulegraph.XXXXXX)
trap 'rm -f "$graph_tmp"' EXIT
PYTHONHASHSEED=0 snakemake all --snakefile analyses/Snakefile \
    --configfile analyses/config/config.yaml "${@}" --rulegraph | dot -Tsvg > "$graph_tmp"
chmod 644 "$graph_tmp"
mv "$graph_tmp" analyses/docs/rulegraph.svg
