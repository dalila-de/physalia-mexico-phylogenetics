#!/usr/bin/env bash
#SBATCH --job-name=physalia_replot
#SBATCH --partition=day
#SBATCH --account=dunn
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=work/ncbi-refinement/provenance/plot-%j.log
set -euo pipefail
export PATH=/vast/palmer/pi/dunn/dd2233/Physalia_Phylogenetics_2026/software/plot_env/bin:/vast/palmer/pi/dunn/dd2233/.conda/envs/IlluminaQC/bin:$PATH
export XDG_CACHE_HOME="$PWD/.cache"
export MPLCONFIGDIR="$PWD/.cache/matplotlib"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
record_dir="analyses/provenance/high-support-node-labels"
mkdir -p "$record_dir"
trap 'rc=$?; echo "$rc" > "$record_dir/exit_code.txt"; exit "$rc"' EXIT
cp work/ncbi-refinement/high_support_plot_batch.sh "$record_dir/command.sh"
snakemake --snakefile analyses/plots.smk --cores 2 --forceall > "$record_dir/workflow.log" 2>&1
