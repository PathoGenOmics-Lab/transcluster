#!/usr/bin/env Rscript

library(tidyverse)
library(logger)

snakemake@source("age_correction.R")

Sys.setlocale("LC_TIME", "English")
theme_set(theme_classic())

agecol <- sym(snakemake@params[["metadata_age_column"]])


log_info("Reading and reformatting haplotype metadata")
metadata <- read_csv(
        snakemake@input[["haplotype_metadata"]],
        col_types = cols(!!agecol := "c")
    ) %>%
    correct.age(., agecol, as.numeric(snakemake@params[["max_age"]])) %>%
    # Tag transmitted samples
    mutate(
        Transmitted = ifelse(
            is.na(cluster_id),
            yes = "No", no = "Yes"
        )
    ) %>%
    # Tag clusters with size >= threshold
    group_by(Haplotype, cluster_id) %>%
    mutate(
        cluster.size.threshold.pass = (n() >= snakemake@params[["cluster_size_threshold"]])
    ) %>%
    ungroup()


log_info("Writing age-corrected haplotype metadata")
metadata %>% write_csv(snakemake@output[["age_corrected_metadata"]])
