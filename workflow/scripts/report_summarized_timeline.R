#!/usr/bin/env Rscript

library(tidyverse)
library(logger)
library(ggpubr)

Sys.setlocale("LC_TIME", "English")
theme_set(theme_pubclean())

datecol <- sym(snakemake@params[["metadata_date_column"]])

log_info("Reading plot data")
plot.data <- lapply(
        snakemake@input[["plot_data_tables"]],
        function(path) read_csv(path, col_types = cols(!!datecol := "D"))
    ) %>%
    bind_rows
n.haplotypes <- length(snakemake@input[["plot_data_tables"]])

log_info("Writing report")
ggplot(plot.data, aes(x = !!datecol, y = Haplotype)) +
    geom_point(
        aes(color = Country),
        size = 1,
        alpha = 0.5,
        position = position_dodge(width = -0.75)
    ) +
    geom_violin(linewidth = 0, fill = "black", alpha = 0.4) +
    scale_x_date(
        date_breaks = "3 month",
        date_minor_breaks = "1 month",
        date_labels = "%b %Y") +
    scale_color_viridis_d()

ggsave(
    snakemake@output[["report"]],
    width = snakemake@params[["width_mm"]],
    height = n.haplotypes * snakemake@params[["height_per_haplotype_mm"]],
    units = "mm",
    limitsize = FALSE
)

log_info("Writing report data")
plot.data %>% write_csv(snakemake@output[["report_data"]])
