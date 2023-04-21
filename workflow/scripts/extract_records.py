#!/usr/bin/env python3

from Bio import SeqIO


def fix_spaces(string):
    return string.replace(" ", snakemake.config["SPACE_REPLACEMENT"])


def parse_record_description(description):
    """If record has GISAID format, parse dates
    returns: record ID, collection date | '', submission date | ''"""
    fields = description.split("|")
    if description.startswith("hCoV-19/") and len(fields) == 3:
        return fields[0], fields[1], fields[2]
    else:
        return fields[0], "", ""


# Read IDs in input file
with open(snakemake.input.ids_file) as f:
    provided_ids = set(line.strip() for line in f)
    n_provided_ids = len(provided_ids)

# Extract records from full FASTA, keeping track of ID modifications
id_equivalence = {}
extracted_records = []
for record in SeqIO.parse(snakemake.params.full_fasta, "fasta"):
    original_record_descr = record.description  # hard copy (str)
    record_id, col_date, sub_date = parse_record_description(original_record_descr)
    if record_id in provided_ids:
        # Replace spaces, remove leading 'hCoV-19/' if present
        if col_date and sub_date:
            new_record_id = fix_spaces("|".join([record_id.replace("hCoV-19/", ""), col_date]))
        else:
            new_record_id = fix_spaces(record_id)
        # Rename ID and keep track of the modification
        record.description = new_record_id
        record.id = new_record_id
        id_equivalence[record_id] = new_record_id
        extracted_records.append(record)

# Write output FASTA
n_records_written = SeqIO.write(
    extracted_records,
    snakemake.output.fasta,
    "fasta"
)
print(f"Written {n_records_written} records out of {n_provided_ids} provided IDs")

# Write table containing IDs/descriptors of extracted records
with open(snakemake.output.extracted_ids, "w") as fw:
    fw.write("original_id,modified_id\n")
    for original_id, modified_id in id_equivalence.items():
        fw.write(f"{original_id},{modified_id}\n")

if n_records_written != n_provided_ids:
    print(f"The following records were not found in '{snakemake.params.full_fasta}':")
    print(*set.symmetric_difference(provided_ids, set(ori for ori,mod in id_equivalence.items())), sep="\n")
