#!/usr/bin/env Rscript

library(tidyverse)
library(logger)
library(ggpubr)

Sys.setlocale("LC_TIME", "English")
theme_set(theme_pubclean())

datecol <- sym(snakemake@params[["metadata_date_column"]])

log_info("Reading metadata")
metadata <- read_csv(
    snakemake@input[["haplotype_metadata"]],
    col_types = cols(!!datecol := "D")
)

log_info("Formatting metadata")
# Extract countries
plot.data <- metadata %>%
    separate(
        Location,
        into = c("Continent", "Country"),
        extra = "drop",
        sep = "\\ /\\ ",
        remove = FALSE
    ) %>%
    drop_na(Haplotype)

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
ggplot(plot.data, aes(x = !!datecol, y = Country)) +
    geom_jitter(
        aes(color = `Top countries`),
        size = 1,
        alpha = 0.75
    ) +
    geom_boxplot(alpha = 0.5) +
    scale_x_date(
        date_breaks = "3 month",
        date_minor_breaks = "1 month",
        date_labels = "%b %Y") +
    scale_color_viridis_d() +
    ggtitle(snakemake@wildcards[["dataset"]])

ggsave(
    snakemake@output[["report"]],
    width = snakemake@params[["width_mm"]],
    height = snakemake@params[["height_mm"]],
    units = "mm"
)

log_info("Writing report data")
plot.data %>% write_csv(snakemake@output[["report_data"]])
