#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO


def format_gisaid_record_id(record_id):
    return record_id.split("|")[0]


Path(snakemake.output.fasta_dir).mkdir(parents=True, exist_ok=True)

ids = set()
for ids_file in Path(snakemake.input.ids_dir).glob("*.txt"):
    haplotype_name = Path(ids_file).stem
    haplotype_record_ids = set()
    with open(ids_file) as f:
        for line in f:
            record_id = line.strip()
            # Add to all ids
            ids.add(record_id)
            # Add to current haplotype ids
            haplotype_record_ids.add(record_id)
    # Write current haplotype ids
    SeqIO.write(
        (record for record in SeqIO.parse(snakemake.params.full_fasta, "fasta") if format_gisaid_record_id(record.id) in haplotype_record_ids),
        Path(snakemake.output.fasta_dir) / f"{haplotype_name}.fasta",
        "fasta"
    )

# Write all ids
records = [record for record in SeqIO.parse(snakemake.params.full_fasta, "fasta") if record.id in ids]
SeqIO.write(
    records,
    snakemake.output.fasta,
    "fasta"
)

print(f"Written {len(records)} out of {len(ids)} records")
