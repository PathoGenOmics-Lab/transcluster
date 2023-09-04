# transcluster

[![PGO badge](https://img.shields.io/badge/PathoGenOmics-Lab-yellow.svg)](https://pathogenomics.github.io/)
[![Snakemake badge](https://img.shields.io/badge/snakemake-≥6.9-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

A Snakemake pipeline to find SARS-CoV-2 clusters in a reference phylogeny.

It requires an input directory (`INPUT_DIR` parameter, default: `input/`) containing at least
one `.txt` file with a list of record IDs (one per line). Each ID list is evaluated
separately, possibly in parallel.
Thus, records listed in a given list of samples are extracted from a FASTA file
containing an aligned sequence dataset (provided through the `ALIGNED_FULL_FASTA` parameter)
and placed on a reference phylogeny in Newick format (`REFERENCE_TREE`parameter ). Then, the
pipeline looks for transmission clusters using a breadth-first search approach, skipping
clades that do not conform to the parametrized clustering criteria (`MIN_SIZE` and `MIN_PROP`
parameters, default: ≥2 members and ≥90% proportion of target samples). Results are written
to a subdirectory within the output directory (`OUTPUT_DIR` parameter, default: `output/`).

All parameters have a default setting except for `ALIGNED_FULL_FASTA`, which must be provided
by the user. Parameters can be set and overridden by modifying the [`config.yaml`](config/config.yaml)
file or through the `--config` command line argument (see the
[Snakemake docs](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#standard-configuration)).

The pipeline can also leverage an HPC environment using SLURM.

## Usage

This pipeline requires `conda`/`mamba` and `snakemake`. The rest of dependencies are
defined in [YAML environment files](workflow/envs) and can be automatically installed
upon execution through the command line argument `--use-conda`.

To run the pipeline with 4 cores:

```bash
snakemake --use-conda -c 4
```

To run a batch job with SLURM (defaults time limit: 15 days):

```bash
sbatch run-slurm.sh
```

Note that the `sbatch` options set in [`run-slurm.sh`](run-slurm.sh) and [`slurm.yaml`](config/slurm.yaml)
might need to be modified to meet the user's needs depending on the particular cluster characteristics.

## Configuration parameters

### `REFERENCE_TREE`

URL or path pointing to a reference phylogeny to place the target samples in.

Select a reference tree with all GISAID sequences up to a date by providing the following `REFERENCE_TREE`
(replace `YYYY`, `MM` and `DD` according to the desired threshold date):

```url
https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/YYYY/MM/DD/public-YYYY-MM-DD.all.masked.pb.gz
```

Alternatively, select a locally stored tree in `.protobuf` format by providing its path.

### `PROBLEMATIC_VCF`

Path or URL pointing to a VCF file containing problematic sites.

### `ALIGNED_FULL_FASTA`

Path for a aligned sequence dataset in FASTA format. For each set of target samples, sequences are sourced
from this file according to their ID. Sequences must be aligned to `REFERENCE_FASTA`.

### `REFERENCE_FASTA`

Path to the reference sequence in FASTA format.

### `INPUT_DIR` and `OUTPUT_DIR`

Path to the input and output directories, respectively.

### `MIN_PROP` and `MIN_SIZE`

Threshold for cluster detections: the minimum proportion of target samples (from 0 to 1, inclusive)
and the minimum number of samples within a cluster, respectively.

### `SPACE_REPLACEMENT`

A string used to replace spaces in FASTA records to avoid truncation of sequence descriptors during
the analysis.

## Testing the pipeline

The pipeline can be tested by running `test-local.sh`.

Cluster search is the most time-consuming task. Run time should be heavily
dependent on tree structure. Time is expected to increase somewhat linearly,
but the diversity of the target set of sequences will be a major factor too.
Expect higher run times with more diverse sets.
