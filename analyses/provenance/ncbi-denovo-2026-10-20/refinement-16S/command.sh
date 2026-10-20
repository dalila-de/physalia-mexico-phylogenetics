#!/usr/bin/env bash
#SBATCH --job-name=physalia_10runs
#SBATCH --partition=day
#SBATCH --account=dunn
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=0-1
#SBATCH --output=work/ncbi-refinement/provenance/slurm-%A_%a.log
set -euo pipefail
markers=(CO1 16S)
marker="${markers[$SLURM_ARRAY_TASK_ID]}"
record_dir="work/ncbi-refinement/provenance/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$record_dir"
trap 'rc=$?; echo "$rc" > "$record_dir/exit_code.txt"; exit "$rc"' EXIT
cp work/ncbi-refinement/batch.sh "$record_dir/command.sh"
sha256sum "analyses/results/alignments/$marker.fasta" > "$record_dir/alignment.sha256"
/vast/palmer/pi/dunn/dd2233/Physalia_Phylogenetics_2026/software/iqtree-3.0.1-Linux/bin/iqtree3 -s "analyses/results/alignments/$marker.fasta" -m MFP -B 1000 --runs 10 -nt 8 -seed 20261012 -pre "work/ncbi-refinement/trees/$marker" > "$record_dir/iqtree.log" 2>&1
