#!/usr/bin/env bash

snakemake \
    --use-conda \
    -c 2 \
    --config ALIGNED_FULL_FASTA="test/full_fasta.fasta" \
             INPUT_DIR="test/input" \
             OUTPUT_DIR="test/output"
