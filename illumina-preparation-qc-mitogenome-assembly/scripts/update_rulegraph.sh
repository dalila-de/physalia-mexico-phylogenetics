#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p docs
graph_tmp=$(mktemp docs/.rulegraph.XXXXXX)
trap 'rm -f "$graph_tmp"' EXIT
PYTHONHASHSEED=0 snakemake all --snakefile Snakefile --configfile config/config.yaml \
    --rulegraph | dot -Tsvg > "$graph_tmp"
mv "$graph_tmp" docs/rulegraph.svg
