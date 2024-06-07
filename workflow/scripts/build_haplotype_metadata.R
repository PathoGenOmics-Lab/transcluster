#!/usr/bin/env Rscript

# Final columns:
# Virus name, Accession ID,
# Patient age, Gender,
# Pango lineage, Collection date, Location,
# Haplotype, cluster_id

library(tidyverse)
library(logger)

METADATA.SELECT.COLS <- c(
    snakemake@params[["metadata_originalid_column"]],
    snakemake@params[["metadata_id_column"]],
    snakemake@params[["metadata_date_column"]],
    snakemake@params[["metadata_location_column"]],
    snakemake@params[["metadata_gender_column"]],
    snakemake@params[["metadata_lineage_column"]],
    snakemake@params[["metadata_age_column"]]
)

log_info("Reading cluster table")
clusters <- read_delim(
    snakemake@input[["clusters"]],
    col_select = c("original_id", "cluster_id")
)

log_info("Reading full metadata and merging")
metadata <- clusters %>%
    left_join(
        read_delim(
            snakemake@input[["full_metadata"]],
            col_select = METADATA.SELECT.COLS
        ),
        by = c("original_id" = snakemake@params[["metadata_originalid_column"]])
    )

log_info("Writing metadata")
metadata %>%
    mutate(Haplotype = snakemake@wildcards[["dataset"]]) %>%
    write_csv(snakemake@output[["metadata"]])