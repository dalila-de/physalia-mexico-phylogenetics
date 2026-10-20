#!/usr/bin/env bash
#SBATCH --job-name=physalia_ncbi
#SBATCH --partition=day
#SBATCH --account=dunn
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=analyses/provenance/slurm-%j.log
set -eo pipefail
module reset
module load miniconda/24.7.1 MAFFT/7.526-GCC-13.3.0-with-extensions
conda activate /vast/palmer/pi/dunn/dd2233/.conda/envs/IlluminaQC
set -u
export PATH=/vast/palmer/pi/dunn/dd2233/Physalia_Phylogenetics_2026/software/iqtree-3.0.1-Linux/bin:$PATH
export XDG_CACHE_HOME="$PWD/.cache"
export MPLCONFIGDIR="$PWD/.cache/matplotlib"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
record_dir="analyses/provenance/slurm-$SLURM_JOB_ID"
mkdir -p "$record_dir"
trap 'rc=$?; echo "$rc" > "$record_dir/exit_code.txt"; exit "$rc"' EXIT
git rev-parse HEAD > "$record_dir/source_commit.txt"
git diff > "$record_dir/source.patch"
tar -czf "$record_dir/source.tar.gz" analyses/Snakefile analyses/workflow analyses/scripts analyses/config/config.yaml references/ncbi/accessions.tsv
conda list --explicit > "$record_dir/environment.txt"
/vast/palmer/pi/dunn/dd2233/Physalia_Phylogenetics_2026/software/plot_env/bin/python -m pip freeze > "$record_dir/analysis_python_packages.txt"
mafft --version > "$record_dir/mafft_version.txt" 2>&1
iqtree3 --version > "$record_dir/iqtree_version.txt"
snakemake all --snakefile analyses/Snakefile --cores 16 --rerun-incomplete --printshellcmds --rerun-triggers mtime --config analysis_python=/vast/palmer/pi/dunn/dd2233/Physalia_Phylogenetics_2026/software/plot_env/bin/python > "$record_dir/workflow.log" 2>&1
