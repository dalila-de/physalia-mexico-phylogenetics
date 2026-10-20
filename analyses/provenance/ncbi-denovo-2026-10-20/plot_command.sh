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
record_dir="analyses/provenance/ncbi-denovo-2026-10-20"
trap 'rc=$?; echo "$rc" > "$record_dir/plot_exit_code.txt"; exit "$rc"' EXIT
cp work/ncbi-refinement/plot_batch.sh "$record_dir/plot_command.sh"
for marker in CO1 16S 18S ITS mitogenome_submitted; do
    python analyses/scripts/report_tree.py "$marker"
done
snakemake --snakefile analyses/plots.smk --cores 2 --forceall > "$record_dir/plot_workflow.log" 2>&1
python analyses/scripts/summarize_run.py > "$record_dir/run_summary.txt"
