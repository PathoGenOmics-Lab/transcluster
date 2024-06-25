#!/usr/bin/env Rscript

library(tidyverse)
library(logger)

Sys.setlocale("LC_TIME", "English")
theme_set(theme_classic())

agecol <- sym(snakemake@params[["metadata_age_column"]])
datecol <- sym(snakemake@params[["metadata_date_column"]])
loccol <- sym(snakemake@params[["metadata_location_column"]])


log_info("Reading age-corrected metadata")
metadata <- lapply(
    snakemake@input[["age_corrected_tables"]],
    function(path) {
        read_csv(
            path,
            col_select = c(Haplotype, !!agecol, !!datecol, !!loccol),
            col_types = cols(!!agecol := "n", !!datecol := "D")
        )
    }
) %>% bind_rows


log_info("Plotting missing ages")
metadata %>%
    mutate(`Missing age` = is.na(!!agecol)) %>%
    ggplot(aes(Haplotype, fill = `Missing age`)) +
        geom_bar(stat = "count") +
        geom_text(
            aes(label = ..count..),
            stat = "count",
            position = position_stack(vjust = 0.5)
        ) +
        scale_fill_viridis_d(begin = 0.2, end = 0.8) +
        ggtitle("Missing age values")

ggsave(
  snakemake@output[["report_age"]],
  width = snakemake@params[["width_mm"]],
  height = snakemake@params[["height_mm"]],
  units = "mm"
)


log_info("Plotting missing locations")
metadata %>%
    mutate(`Missing location` = is.na(!!loccol)) %>%
    ggplot(aes(Haplotype, fill = `Missing location`)) +
        geom_bar(stat = "count") +
        geom_text(
            aes(label = ..count..),
            stat = "count",
            position = position_stack(vjust = 0.5)
        ) +
        scale_fill_viridis_d(begin = 0.2, end = 0.8) +
        ggtitle("Missing location values")

ggsave(
  snakemake@output[["report_location"]],
  width = snakemake@params[["width_mm"]],
  height = snakemake@params[["height_mm"]],
  units = "mm"
)

log_info("Plotting missing dates")
metadata %>%
    mutate(`Missing date` = is.na(!!datecol)) %>%
    ggplot(aes(Haplotype, fill = `Missing date`)) +
        geom_bar(stat = "count") +
        geom_text(
            aes(label = ..count..),
            stat = "count",
            position = position_stack(vjust = 0.5)
        ) +
        scale_fill_viridis_d(begin = 0.2, end = 0.8) +
        ggtitle("Missing date values")

ggsave(
  snakemake@output[["report_date"]],
  width = snakemake@params[["width_mm"]],
  height = snakemake@params[["height_mm"]],
  units = "mm"
)
