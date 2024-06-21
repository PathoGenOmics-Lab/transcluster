#!/usr/bin/env bash

library(tidyverse)
library(logger)
library(ggpubr)

snakemake@source("age_correction.R")

Sys.setlocale("LC_TIME", "English")
theme_set(theme_classic())

idcol <- sym(snakemake@params[["metadata_id_column"]])
agecol <- sym(snakemake@params[["metadata_age_column"]])
datecol <- sym(snakemake@params[["metadata_date_column"]])

log_info("Reading age-corrected haplotype IDs, ages and dates")
age.metadata <- read_csv(
    snakemake@input[["age_corrected_metadata"]],
    col_select = c(!!idcol, !!agecol, !!datecol),
    col_types = cols(!!agecol := "n", !!datecol := "D")
)

min.haplotype.date <- min(
    age.metadata[[snakemake@params[["metadata_date_column"]]]], na.rm = TRUE
)
max.haplotype.date <- max(
    age.metadata[[snakemake@params[["metadata_date_column"]]]], na.rm = TRUE
)

log_info("Reading background IDs, ages and dates, filtering by date, and correcting ages")
bg.metadata <- read_delim(
        snakemake@input[["full_metadata"]],
        col_select = c(!!idcol, !!agecol, !!datecol),
        col_types = cols(!!agecol := "c", !!datecol := "D")
    ) %>%
    filter(min.haplotype.date <= !!datecol, !!datecol <= max.haplotype.date) %>%
    correct.age(., agecol)

log_info("Building combined age table")
# The pipeline ensures age.metadata is a subset of bg.metadata
plot.data <- bg.metadata %>%
    mutate(
        is.haplotype = !!idcol %in% age.metadata[[idcol]]
    )

log_info("Building age report")
ks <- ks.test(
    plot.data %>% filter(is.haplotype) %>% pull(!!agecol),
    plot.data %>% filter(!is.haplotype) %>% pull(!!agecol)
)
plot.data %>%
    ggplot(aes(is.haplotype, !!agecol)) +
    stat_compare_means(method = "wilcox.test", paired = FALSE) +
    ggtitle(
        snakemake@wildcards[["dataset"]],
        paste0(ks$method, ": D = ", ks$statistic, "; p = ", ks$p.value)
    )
ggsave(
    snakemake@output[["report"]],
    width = snakemake@params[["width_mm"]],
    height = snakemake@params[["height_mm"]],
    units = "mm"
)

log_info("Writing plot data")
plot.data %>% write_csv(snakemake@output[["report_data"]])
