#!/usr/bin/env Rscript

library(tidyverse)
library(logger)
library(ggpubr)

Sys.setlocale("LC_TIME", "English")
theme_set(theme_pubclean())

datecol <- data_sym(snakemake@params[["metadata_date_column"]])

log_info("Reading metadata")
metadata <- read_csv(snakemake@input[["haplotype_metadata"]])

log_info("Formatting metadata")
plot.data <- metadata %>%
    separate(
        Location,
        into = c("Continent", "Country"),
        extra = "drop",
        sep = "\\ /\\ ",
        remove = FALSE
    ) %>%
    drop_na(Haplotype) %>%
    # Rank countries
    group_by(Country) %>%
    mutate(
        Country.Rank = rank(n(), ties.method = "first"),
        Country = ifelse(Country.Rank <= 10, Country, "[Other]")
    ) %>%
    ungroup()

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
    height = snakemake@params[["height_mm"]],
    units = "mm"
)

log_info("Writing report data")
plot.data %>% write_csv(snakemake@output[["report_data"]])
