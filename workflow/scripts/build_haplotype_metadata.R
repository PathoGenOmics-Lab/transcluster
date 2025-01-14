#!/usr/bin/env Rscript

# Input clusters columns: cluster_id, label

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

log_info("Reading cluster table and original IDs")
clustering.results <- read_delim(snakemake@input[["clusters"]])
if (is_empty(clustering.results)) {
    clustering.results <- tibble(
        label = character(),
        cluster_id = character()
    )
}
clusters <- right_join(
    clustering.results,
    read_delim(snakemake@input[["extracted_ids"]]),
    by = c("label" = "modified_id")
)

log_info("Reading full metadata and merging")
metadata <- clusters %>%
    left_join(
        read_delim(
            snakemake@input[["full_metadata"]],
            col_select = all_of(METADATA.SELECT.COLS)
        ),
        by = c("original_id" = snakemake@params[["metadata_originalid_column"]])
    )
# Ensure all columns are present
metadata <- metadata[, METADATA.SELECT.COLS]

log_info("Writing metadata")
metadata %>%
    mutate(Haplotype = snakemake@wildcards[["dataset"]]) %>%
    write_csv(snakemake@output[["metadata"]])
