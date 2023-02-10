#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO


def format_gisaid_record_id(record_id):
    return record_id.split("|")[0]


haplotype_name = Path(snakemake.input.ids_dir).stem
ids = set()
with open(ids_file) as f:
    for line in f:
        record_id = line.strip()
        # Add to all ids
        ids.add(record_id)

# Write all ids
SeqIO.write(
    (record for record in SeqIO.parse(snakemake.params.full_fasta, "fasta") if record.id in ids),
    snakemake.output.fasta,
    "fasta"
)
