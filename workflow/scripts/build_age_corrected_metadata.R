#!/usr/bin/env Rscript

library(tidyverse)
library(logger)

Sys.setlocale("LC_TIME", "English")
theme_set(theme_classic())

agecol <- sym(snakemake@params[["metadata_age_column"]])


log_info("Reading and reformatting haplotype metadata")
metadata <- read_csv(
        snakemake@input[["haplotype_metadata"]],
        col_types = cols(!!agecol := "c")
    ) %>%
    # Try to fix ages
    mutate(
        !!agecol := case_when(
        # e.g. 50
        !is.na(as.numeric(!!agecol)) ~
            as.numeric(!!agecol),
        # e.g. 2 months / 1 month / 6 mos
        grepl("^[0-9]+ ([mM]onths?|mos?)$", !!agecol) ~
            parse_number(!!agecol) / 12,
        # e.g. 1 year / 2 Years
        grepl("^[0-9]+ [yY]ears?$", !!agecol) ~
            parse_number(!!agecol),
        # e.g. 2 years 3 months
        grepl("^[0-9]+ [yY]ears? [0-9]+ ([mM]onths?|mos?)$", !!agecol) ~
            as.numeric(
                gsub("^([0-9]+) [yY]ears? [0-9]+ ([mM]onths?|mos?)$", "\\1", !!agecol)
            ) +
            as.numeric(
                gsub("^[0-9]+ [yY]ears? ([0-9]+) ([mM]onths?|mos?)$", "\\1", !!agecol)
            ) / 12,
        ),
        # all the rest
        TRUE ~ NA
    ) %>%
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
