#!/usr/bin/env bash
#SBATCH -J tc
#SBATCH --part global
#SBATCH --qos medium
#SBATCH -t 7-00:00:00
#SBATCH --mem 4GB
#SBATCH -c 1
#SBATCH -o slurm-transcluster-%j_main.out

SHADOW="/scr/$USER"

snakemake \
    --use-conda \
    --shadow-prefix "$SHADOW" \
    --workflow-profile profiles/garnatxa
