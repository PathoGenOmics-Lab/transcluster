#!/usr/bin/env bash
#SBATCH -J tc
#SBATCH --part global
#SBATCH --qos medium
#SBATCH -t 7-00:00:00
#SBATCH --mem 4GB
#SBATCH -c 1
#SBATCH -o slurm-transcluster-%j_main.out

MAXJOBS=4
SLURM_CONFIG="sm/slurm.yaml"

snakemake \
    --use-conda \
    --cluster-config "$SLURM_CONFIG" \
    --cluster "sbatch -J $SLURM_JOB_NAME-{cluster.name} -t {cluster.time} --qos {cluster.qos} -N {cluster.nodes} --cpus-per-task {cluster.cpus} --mem {cluster.mem} --output {cluster.out}" \
    --jobs $MAXJOBS
