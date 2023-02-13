
"""
input dir:  input
output dir: output/clusters
"""

import re
import json
import requests
from pathlib import Path

configfile: "sm/config.yaml"


def download_file(url, path):
    response = requests.get(url)
    with open(path, "wb") as fw:
        fw.write(response.content)


DATASETS = [path.stem for path in Path("input").glob("*.txt")]


rule all:
    input:
        expand("output/{dataset}/clusters", dataset=DATASETS)


rule download_problematic_vcf:
    threads: 1
    output:
        vcf_path = "data/problematic_sites_sarsCov2.vcf"
    run:
        download_file(config["PROBLEMATIC_VCF_URL"], output.vcf_path)


rule download_reference_tree:
    threads: 1
    output:
        tree_path = "data/public.all.masked.pb.gz"
    run:
        download_file(config["REFERENCE_TREE_URL"], output.tree_path)


rule extract_records:
    threads: 1
    conda: "sm/bio_env.yaml"
    input:
        ids_file = "input/{dataset}.txt"
    params:
        full_fasta = config["ALIGNED_FULL_FASTA"]
    output:
        fasta = "output/{dataset}/sequences.aligned.fasta"
    script:
        "scripts/extract_records.py"


rule phylogenetic_placement:
    threads: 32
    shadow: "shallow"
    conda: "sm/usher_env.yaml"
    input:
        new_samples = "output/{dataset}/sequences.aligned.fasta",
        masking_vcf = "data/problematic_sites_sarsCov2.vcf",
        tree_protobuf = "data/public.all.masked.pb.gz"
    params:
        reference_fasta = "data/reference.fasta",
        reference_id = "NC_045512.2",
        output_directory = "pp_out"  # if shadow == "shallow", it will be removed
    output:
        tree = "output/{dataset}/tree.nwk"
    shell:
        """
        cat {params.reference_fasta} {input.new_samples} > input.fasta
        faToVcf -maskSites={input.masking_vcf} -ref={params.reference_id} input.fasta input.vcf
        usher -T {threads} -i {input.tree_protobuf} -v input.vcf -u -d {params.output_directory}
        cp {params.output_directory}/uncondensed-final-tree.nh {output.tree}
        """


rule build_phylo4:
    threads: 1
    conda: "sm/r_env.yaml"
    input:
        tree = "output/{dataset}/tree.nwk"
    output:
        tree_p4 = "output/{dataset}/phylo4.RData"
    script:
        "scripts/build_phylo4.R"


rule calculate_clusters:
    threads: 1
    conda: "sm/r_env.yaml"
    input:
        ids_file = "input/{dataset}.txt",
        tree_p4 = "output/{dataset}/phylo4.RData"
    params:
        min_prop = 0.9,
        min_size = 2
    output:
        out_dir = directory("output/{dataset}/clusters"),
    script:
        "scripts/phylo_clusters.R"
