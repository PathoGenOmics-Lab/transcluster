#!/usr/bin/env Rscript

library(tidyverse)
library(logger)
library(ggpubr)

Sys.setlocale("LC_TIME", "English")
theme_set(theme_pubclean())

datecol <- sym(snakemake@params[["metadata_date_column"]])
loccol <- sym(snakemake@params[["metadata_location_column"]])

log_info("Reading plot data")
plot.data <- lapply(
        snakemake@input[["haplotype_metadata_tables"]],
        function(path) read_csv(path, col_types = cols(!!datecol := "D"))
    ) %>%
    bind_rows
n.haplotypes <- length(snakemake@input[["haplotype_metadata_tables"]])

log_info("Formatting metadata")
if (snakemake@params$separate_location) {
    # Extract countries
    plot.data <- plot.data %>%
        separate(
            !!loccol,
            into = c("Continent", "Country"),
            extra = "drop",
            sep = "\\ /\\ ",
            remove = FALSE
        ) %>%
        drop_na(Haplotype)
} else {
    # Rename country column
    plot.data <- plot.data %>%
        rename(Country = !!loccol)
        drop_na(Haplotype)
}

# Rank countries
top.countries <- plot.data %>%
  count(Country, sort = TRUE) %>%
  top_n(snakemake@params[["n_top_countries"]], n) %>%
  pull(Country)

# Add column for top countries
plot.data <- plot.data %>%
    mutate(
        `Top countries` = ifelse(
            Country %in% top.countries,
            Country, "[Other]"
        )
    )

log_info("Writing report")
ggplot(plot.data, aes(x = !!datecol, y = Haplotype)) +
    geom_violin(linewidth = 0, fill = "black", alpha = 0.4) +
    geom_point(
        aes(color = `Top countries`),
        size = 2,
        alpha = 0.25,
        position = position_dodge(width = -0.75)
    ) +
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
