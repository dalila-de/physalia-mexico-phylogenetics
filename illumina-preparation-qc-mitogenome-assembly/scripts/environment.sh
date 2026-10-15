#!/usr/bin/env bash
# Source from the repository root on McCleary.
module reset
module load miniconda/24.7.1 Trimmomatic/0.39-Java-11 FastQC/0.12.1-Java-11 \
    Kraken2/2.1.3-gompi-2020b BWA/0.7.17-GCC-10.2.0 SAMtools/1.16.1-GCCcore-10.2.0
conda activate /vast/palmer/pi/dunn/dd2233/.conda/envs/IlluminaQC
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$PWD/.cache}"
