#!/usr/bin/env bash
#SBATCH --job-name=physalia_illumina_qc
#SBATCH --partition=ycga
#SBATCH --account=dunn
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00

set -eo pipefail
source scripts/environment.sh
set -u
repo_dir="$PWD"
run_dir=$(python -c 'import yaml; from pathlib import Path; c=yaml.safe_load(open("config/config.yaml")); print(Path(c["scratch"])/c["run_id"])')
record_dir="$repo_dir/provenance/slurm-${SLURM_JOB_ID}"
mkdir -p "$run_dir" "$record_dir"
export XDG_CACHE_HOME="$run_dir/.cache"
export GETORG_PATH=$(python -c 'import yaml; print(yaml.safe_load(open("config/config.yaml"))["getorganelle_db"])')
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
trap 'rc=$?; echo "$rc" > "$record_dir/exit_code.txt"; python "$repo_dir/scripts/status.py" > "$record_dir/final_status.tsv"; exit "$rc"' EXIT

git diff --exit-code
git diff --cached --exit-code
git rev-parse HEAD > "$record_dir/source_commit.txt"
git archive HEAD | gzip > "$record_dir/source.tar.gz"
cp config/config.yaml "$record_dir/config.yaml"
module -t list > "$record_dir/modules.txt" 2>&1
conda list --explicit > "$record_dir/conda-explicit.txt"
{
    python --version
    snakemake --version
    get_organelle_from_reads.py --version
    kraken2 --version
    samtools --version
    fastqc --version
} > "$record_dir/tool_versions.txt" 2>&1
python scripts/preflight.py --output "$record_dir/preflight.json"
get_organelle_config.py --check --use-version 0.0.1 --config-dir "$GETORG_PATH" > "$record_dir/getorganelle_db_check.txt" 2>&1
python scripts/record_references.py --output "$record_dir/references.tsv"

run_snakemake() {
    snakemake "$@" --snakefile "$repo_dir/Snakefile" --directory "$run_dir" \
        --configfile "$repo_dir/config/config.yaml" --executor local \
        --cores "$SLURM_CPUS_PER_TASK" --resources mem_mb=250000 \
        --rerun-incomplete --keep-going --scheduler greedy --printshellcmds
}

# A failed pilot stops here. Successful pilot products also count toward all.
run_snakemake pilot > "$record_dir/pilot.log" 2>&1
python scripts/status.py --require-pilot > "$record_dir/pilot_status.tsv"

# Bound retained temporary reads by processing at most four specimens at a time.
mapfile -t targets < <(python - <<'PY'
import csv
from pathlib import Path
import yaml
c = yaml.safe_load(open('config/config.yaml'))
for row in csv.DictReader(open(c['samples']), delimiter='\t'):
    if row['sample'] not in c['pilot']:
        print(Path.cwd() / 'results' / c['run_id'] / row['sample'] / 'summary.json')
PY
)
failed=0
for ((offset=0; offset<${#targets[@]}; offset+=4)); do
    if ! run_snakemake "${targets[@]:offset:4}" > "$record_dir/batch-${offset}.log" 2>&1; then
        failed=1
    fi
done
if ((failed)); then
    exit 1
fi
run_snakemake all > "$record_dir/cohort.log" 2>&1
