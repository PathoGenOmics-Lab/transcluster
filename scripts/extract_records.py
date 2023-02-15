#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO


def remove_dates_gisaid_record_id(record_id):
    return record_id.split("|")[0]


def iter_format_records(records):
    for record in records:
        record_id_no_dates = remove_dates_gisaid_record_id(record.description)
        if record_id_no_dates in ids:
            fixed_record_id = record_id_no_dates.replace(" ", snakemake.params.space_replacement)
            record.description = fixed_record_id
            record.id = fixed_record_id
            yield record


ids = set()
with open(snakemake.input.ids_file) as f:
    for line in f:
        record_id = line.strip()
        # Add to all ids
        ids.add(record_id)


# Write all ids
n_records_written = SeqIO.write(
    iter_format_records(SeqIO.parse(snakemake.params.full_fasta, "fasta")),
    snakemake.output.fasta,
    "fasta"
)
print(f"Written {n_records_written} records out of {len(ids)}")
