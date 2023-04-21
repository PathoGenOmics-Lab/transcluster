#!/usr/bin/env bash
#SBATCH -J tc
#SBATCH --part global
#SBATCH --qos long
#SBATCH -t 15-00:00:00
#SBATCH --mem 4GB
#SBATCH -c 1
#SBATCH -o slurm-transcluster-%j_main.out

MAXJOBS=4
SHADOW="/scr/$USER"
SLURM_CONFIG="config/slurm.yaml"

snakemake \
    --use-conda \
    --shadow-prefix "$SHADOW" \
    --cluster-config "$SLURM_CONFIG" \
    --cluster "sbatch -J $SLURM_JOB_NAME-{cluster.name} -t {cluster.time} --qos {cluster.qos} -N {cluster.nodes} --cpus-per-task {cluster.cpus} --mem {cluster.mem} --output {cluster.out}" \
    --jobs $MAXJOBS
