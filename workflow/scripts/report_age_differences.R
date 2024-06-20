#!/usr/bin/env Rscript

library(tidyverse)
library(logger)

Sys.setlocale("LC_TIME", "English")
theme_set(theme_classic())

combinations <- function(elements) {
  # Generate all possible pairs of elements
  if (length(elements) < 2) {
    return(list(elements))
  }
  pairs <- combn(elements, 2)
  combos <- list()
  for (i in seq_len(ncol(pairs))) {
    combos[[i]] <- pairs[, i]
  }
  return(combos)
}

agecol <- sym(snakemake@params[["metadata_age_column"]])


log_info("Reading plot data and fixing ages")
age.metadata <- lapply(
        snakemake@input[["haplotype_metadata_tables"]],
        function(path) read_csv(path, col_types = cols(!!agecol := "c"))
    ) %>%
    bind_rows %>%
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


log_info("Plotting age differences between haplotypes")
age.metadata %>%
    ggplot(aes(x = Haplotype, y = !!agecol)) +
    geom_jitter(alpha = 0.5) +
    geom_violin() +
    geom_boxplot(alpha = 0.2, varwidth = TRUE) +
    stat_summary(fun.y = mean, geom = "point") +
    stat_compare_means(
        comparisons = combinations(unique(age.metadata$Haplotype)),
        method = "wilcox", paired = FALSE
    )

ggsave(
    snakemake@output[["report_haplotype"]],
    width = snakemake@params[["width_mm"]],
    height = snakemake@params[["height_mm"]],
    units = "mm"
)


log_info("Plotting age differences regarding transmission")
age.metadata %>%
    filter(cluster.size.threshold.pass) %>%
    ggplot(aes(x = Haplotype, y = !!agecol, fill = Transmitted)) +
    geom_jitter(alpha = 0.5) +
    geom_violin() +
    geom_boxplot(alpha = 0.2, varwidth = TRUE) +
    stat_summary(fun.y = mean, geom = "point") +
    stat_compare_means(method = "wilcox", paired = FALSE) +
    scale_fill_viridis_d()

ggsave(
    snakemake@output[["report_transmission"]],
    width = snakemake@params[["width_mm"]],
    height = snakemake@params[["height_mm"]],
    units = "mm"
)


log_info("Writing age-corrected haplotype metadata")
age.metadata %>% write_csv(snakemake@output[["report_data"]])
