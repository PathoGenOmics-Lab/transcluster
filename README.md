# transcluster

A Snakemake pipeline to find SARS-CoV-2 clusters in a reference phylogeny.


## Usage

Command line (with 4 cores):

```bash
snakemake --use-conda -c 4
```

You can modify the default configuration in `sm/config.yaml` like this:

```bash
snakemake --use-conda -c 4 --config ALIGNED_FULL_FASTA="my_full.aligned.fasta" OUTPUT_DIR="my_output_dir"
```

Run defalult batch job with slurm (7 days tops):

```bash
sbatch run-slurm.sh
```


## Some configuration details

### `REFERENCE_TREE`

You can select a reference tree with all GISAID sequences up to a date (e.g. yyyy/mm/dd) by providing the following `REFERENCE_TREE_URL`:

```
https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/{yyyy}/{mm}/{dd}/public-{yyyy}-{mm}-{dd}.all.masked.pb.gz
```

You can also select a locally stored protobuf tree by providing its path.

### `ALIGNED_FULL_FASTA`

It must be aligned to `REFERENCE_FASTA`!

### `PROBLEMATIC_VCF`

You probably don't need to change it, but you may want to be more or less flexible regarding problematic sites.
To do that, you can select a locally stored VCF with problematic sites by providing its path.
