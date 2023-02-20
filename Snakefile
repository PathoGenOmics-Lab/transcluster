
"""
input dir:  input
output dir: output/clusters
"""

import re
import json
import requests
from pathlib import Path

configfile: "sm/config.yaml"


def read_fasta_IDs(path):
    ids = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split(" ")[0])
    return ids


def download_file(url, path):
    response = requests.get(url)
    with open(path, "wb") as fw:
        fw.write(response.content)


def copy_file(from_path, to_path, chunk_size=1024):
    with open(from_path, "rb") as f:
        with open(to_path, "wb") as fw:
            fw.write(f.read())


INPUT_DIR = Path(config["INPUT_DIR"])
OUTPUT_DIR = Path(config["OUTPUT_DIR"])
DATASETS = [path.stem for path in INPUT_DIR.glob("*.txt")]


rule all:
    input:
        expand(OUTPUT_DIR/"{dataset}/clusters", dataset=DATASETS)


rule get_problematic_vcf:
    threads: 1
    shadow: "shallow"
    output:
        vcf_path = "data/problematic_sites_sarsCov2.vcf"
    run:
        if config["PROBLEMATIC_VCF"].startswith("https://"):
            download_file(config["PROBLEMATIC_VCF"], output.vcf_path)
        else:
            copy_file(config["PROBLEMATIC_VCF"], output.vcf_path)


rule get_reference_tree:
    threads: 1
    shadow: "shallow"
    output:
        tree_path = "data/public.all.masked.pb.gz"
    run:
        if config["REFERENCE_TREE"].startswith("https://"):
            download_file(config["REFERENCE_TREE"], output.tree_path)
        else:
            copy_file(config["REFERENCE_TREE"], output.tree_path)


rule extract_records:
    threads: 1
    shadow: "shallow"
    conda: "sm/bio_env.yaml"
    input:
        ids_file = INPUT_DIR/"{dataset}.txt"
    params:
        full_fasta = config["ALIGNED_FULL_FASTA"]
    output:
        fasta = OUTPUT_DIR/"{dataset}/sequences.aligned.fasta",
        extracted_ids = OUTPUT_DIR/"{dataset}/extracted_ids.csv"
    script:
        "scripts/extract_records.py"


rule phylogenetic_placement:
    threads: 32
    shadow: "shallow"
    conda: "sm/usher_env.yaml"
    input:
        new_samples = OUTPUT_DIR/"{dataset}/sequences.aligned.fasta",
        masking_vcf = "data/problematic_sites_sarsCov2.vcf",
        tree_protobuf = "data/public.all.masked.pb.gz"
    params:
        reference_fasta = config["REFERENCE_FASTA"],
        reference_id = read_fasta_IDs(config["REFERENCE_FASTA"])[0],
        output_directory = "pp_out"  # if shadow == "shallow", it will be removed
    output:
        tree = OUTPUT_DIR/"{dataset}/tree.nwk"
    shell:
        """
        cat {params.reference_fasta} {input.new_samples} > input.fasta
        faToVcf -maskSites='{input.masking_vcf}' -ref='{params.reference_id}' input.fasta input.vcf
        usher -T {threads} -i {input.tree_protobuf} -v input.vcf -u -d {params.output_directory}
        cp {params.output_directory}/uncondensed-final-tree.nh {output.tree}
        """


rule build_phylo4:
    threads: 1
    shadow: "shallow"
    conda: "sm/r_env.yaml"
    input:
        tree = OUTPUT_DIR/"{dataset}/tree.nwk",
        extracted_ids = OUTPUT_DIR/"{dataset}/extracted_ids.csv"
    output:
        tree_p4 = OUTPUT_DIR/"{dataset}/phylo4.RData"
    script:
        "scripts/build_phylo4.R"


rule calculate_clusters:
    threads: 1
    shadow: "shallow"
    conda: "sm/r_env.yaml"
    input:
        tree_p4 = OUTPUT_DIR/"{dataset}/phylo4.RData",
        extracted_ids = OUTPUT_DIR/"{dataset}/extracted_ids.csv"
    params:
        min_prop = 0.9,
        min_size = 2,
        log_every_seconds = 60
    output:
        out_dir = directory(OUTPUT_DIR/"{dataset}/clusters")
    script:
        "scripts/phylo_clusters.R"
