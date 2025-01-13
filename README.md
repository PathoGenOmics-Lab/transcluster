# transcluster

[![PGO badge](https://img.shields.io/badge/PathoGenOmics-Lab-yellow.svg)](https://pathogenomics.github.io/)
[![Release](https://img.shields.io/github/release/PathoGenOmics-Lab/transcluster.svg)](https://github.com/PathoGenOmics-Lab/transcluster/releases)
[![Snakemake badge](https://img.shields.io/badge/snakemake-≥8.12-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

A Snakemake pipeline to find SARS-CoV-2 clusters in a reference phylogeny.

![Simplified workflow rules](/docs/rules-simple.png)

It requires an input directory (`INPUT_DIR` parameter, default: `input/`) containing at least
one `.txt` file with a list of sample IDs (one per line). Each ID list is evaluated
separately, possibly in parallel.

The sequence records from each list of target samples are extracted from a FASTA file
containing an aligned sequence dataset (provided through the `ALIGNED_FULL_FASTA` parameter),
and placed on a reference phylogeny in Newick format (`REFERENCE_TREE` parameter), using [UShER](https://usher-wiki.readthedocs.io/en/latest/UShER.html).
Sample identifiers are transformed beforehand so that samples are deduplicated automatically.

Then, the
pipeline looks for clusters using a breadth-first search approach, skipping
clades that do not conform to the parametrized clustering criteria (`MIN_SIZE` and `MIN_PROP`
parameters, default: ≥2 members and ≥90% proportion of target samples),
using [`tcfinder`](https://github.com/PathoGenOmics-Lab/tcfinder).

Additionally,
a fitness metric is estimated for each cluster by calculating the the ratio of the number
of tips in the cluster to the number of sequences deposited in GISAID during the same time period.
These "background" samples include those collected between the first and last case within the
cluster. The time window is symmetrically padded with a specified number of days (`FITNESS_PADDING_DAYS`
parameter, default: 0 and 7). The denominator helps control for the effect of uneven sequencing
efforts across different times and locations, and is calculated in two ways:

1. *Adding* the number of background samples from all the countries where the cluster was observed.
2. *Slicing* the cluster into each country where the cluster was observed, separately calculating the background
   samples for each country with its different time windows, and then adding them together.

![Estimated fitness denominator](/docs/estimated-fitness-denominator.png)

Finally, a series of figures (including sample timelines, comparisons of estimated fitness and missing data),
are generated. Results are written to a subdirectory within the output directory
(`OUTPUT_DIR` parameter, default: `output/`).

All parameters have a default value except for `ALIGNED_FULL_FASTA` and `METADATA`, which must
be provided by the user. Parameters can be set and overridden by modifying the
[`config.yaml`](config/config.yaml) file or through the `--config` command line argument (see
the [Snakemake docs](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#standard-configuration)).

## Usage

### Local execution

The workflow requires `conda`/`mamba` and `snakemake`. The rest of dependencies are
defined in [YAML environment files](workflow/envs) and can be automatically installed
upon execution through the command line argument `--use-conda`.

To run the pipeline using 4 cores, with the [`default` workflow profile](/profiles/default):

```bash
snakemake --use-conda -c 4
```

### Cluster or cloud execution

The workflow can be easily run on a cluster or in the cloud through
[Snakemake executor plugins](https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#cluster-or-cloud-execution).
Further configuration can be set using [Snakemake profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).

We provide the [`run-slurm.sh`](/run-slurm.sh) Slurm batch script as an example.
It requires the [`slurm` executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).
Note that the `SBATCH` options and the profile definition might need to be
modified to meet the user's needs depending on the particular cluster.

## Configuration parameters

### `ALIGNED_FULL_FASTA`

Path to the aligned sequence dataset in FASTA format. For each set of target samples, sequences are sourced
from this file according to their ID. Sequences must be aligned to `REFERENCE_FASTA`.

### `METADATA`

Path to the GISAID metadata of the input sequences in tabular format.

### `METADATA_COLS / FULL_METADATA`

This section defines the GISAID metadata column names. GISAID may update their column names over time,
so these settings ensure compatibility.

- `ORIGINAL_ID`: the sample identifier as it appears in the sequence record in the FASTA file (default: "Virus name").
- `UNIQUE_ID`: a unique identifier guaranteed to be specific to each sample (default: "Accession ID", the EPI_ISL IDs).
- `DATE`: date when the sample was collected (default: "Collection date").
- `LOCATION`: location where the sample was collected (default: "Location").
- `GENDER`: gender of the host or patient from whom the sample was collected (default: "Gender").
- `LINEAGE`: viral lineage assigned to the sample (default: "Pango lineage").
- `AGE`: age of the host or patient from whom the sample was collected (default: "Patient age").

### `REFERENCE_TREE`

URL or path pointing to a reference phylogeny to place the target samples in.

Select a reference tree containing the appropriate GISAID sequences up to a date
by providing the following `REFERENCE_TREE` (replace `YYYY`, `MM` and `DD`
according to the desired threshold date):

```url
https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/YYYY/MM/DD/public-YYYY-MM-DD.all.masked.pb.gz
```

Alternatively, select a locally stored tree in
[UShER](https://usher-wiki.readthedocs.io/en/latest/UShER.html) protobuf format by providing its path.

### `PROBLEMATIC_VCF`

Path or URL pointing to a VCF file containing problematic sites.

### `REFERENCE_FASTA`

Path to the reference sequence in FASTA format.

### `INPUT_DIR` and `OUTPUT_DIR`

Path to the input and output directories, respectively.

### `MIN_PROP` and `MIN_SIZE`

Thresholds for cluster detections: the minimum proportion of target samples (from 0 to 1, inclusive)
and the minimum number of samples within a cluster, respectively.

### `SPACE_REPLACEMENT`

A string used to replace spaces in FASTA records to avoid truncation of sequence descriptors during
the analysis.

### `FITNESS_PADDING_DAYS`

A list of integers that specify the number of days to pad the fitness estimation denominator time windows with.

## Testing the pipeline

The workflow can be run with test data by executing [`test-local.sh`](/test-local.sh).
