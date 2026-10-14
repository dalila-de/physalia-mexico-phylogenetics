#!/usr/bin/env bash
#SBATCH --job-name=physalia_sanger
#SBATCH --partition=day
#SBATCH --account=dunn
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --chdir=/vast/palmer/pi/dunn/dd2233/physalia-mexico-phylogenetics
#SBATCH --output=/vast/palmer/pi/dunn/dd2233/physalia-mexico-phylogenetics/sanger-2026-10-14/provenance/slurm-%j.log
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
record_dir="sanger-2026-10-14/provenance/slurm-$SLURM_JOB_ID"
mkdir -p "$record_dir"
trap 'rc=$?; echo "$rc" > "$record_dir/exit_code.txt"; exit "$rc"' EXIT
git rev-parse HEAD > "$record_dir/source_commit.txt"
git diff > "$record_dir/source.patch"
cp sanger-2026-10-14/scripts/*.py sanger-2026-10-14/Snakefile "$record_dir/"
/vast/palmer/pi/dunn/dd2233/Physalia_Phylogenetics_2026/software/plot_env/bin/python -m pip freeze > "$record_dir/analysis_python_packages.txt"
conda list --explicit > "$record_dir/environment.txt"
mafft --version > "$record_dir/mafft_version.txt" 2>&1
iqtree3 --version > "$record_dir/iqtree_version.txt"
snakemake all --snakefile sanger-2026-10-14/Snakefile --cores 8 --config analysis_python=/vast/palmer/pi/dunn/dd2233/Physalia_Phylogenetics_2026/software/plot_env/bin/python --rerun-incomplete --printshellcmds > "$record_dir/workflow.log" 2>&1
