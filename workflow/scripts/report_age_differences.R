#!/usr/bin/env Rscript

library(tidyverse)
library(logger)
library(ggpubr)

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


log_info("Reading age-corrected metadata")
age.metadata <- lapply(
  snakemake@input[["age_corrected_tables"]],
  function(path) {
    read_csv(
      path,
      col_select = c(
        Haplotype, !!agecol, Transmitted, cluster.size.threshold.pass
      ),
      col_types = cols(!!agecol := "n")
    )
  }
) %>% bind_rows
n.haplotypes <- length(snakemake@input[["age_corrected_tables"]])

log_info("Plotting age differences between haplotypes")
p <- age.metadata %>%
    ggplot(aes(x = Haplotype, y = !!agecol)) +
    geom_jitter(alpha = 0.2) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(alpha = 0.5) +
    stat_summary(fun.y = mean, geom = "point")
p

ggsave(
    snakemake@output[["report_haplotype_simple"]],
    width = n.haplotypes * snakemake@params[["width_per_haplotype_mm"]],
    height = snakemake@params[["height_mm"]],
    units = "mm"
)

p +
    stat_compare_means(
        comparisons = combinations(unique(age.metadata$Haplotype)),
        method = "wilcox", paired = FALSE
    )

ggsave(
    snakemake@output[["report_haplotype"]],
    width = n.haplotypes * snakemake@params[["width_per_haplotype_mm"]],
    height = snakemake@params[["height_mm"]],
    units = "mm"
)


log_info("Plotting age differences regarding transmission")
age.metadata %>%
    filter(cluster.size.threshold.pass) %>%
    ggplot(aes(x = Haplotype, y = !!agecol, fill = Transmitted)) +
    geom_jitter(alpha = 0.2) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(alpha = 0.5) +
    stat_summary(fun.y = mean, geom = "point") +
    stat_compare_means(method = "wilcox", paired = FALSE) +
    scale_fill_viridis_d()

ggsave(
    snakemake@output[["report_transmission"]],
    width = n.haplotypes * snakemake@params[["width_per_haplotype_mm"]],
    height = snakemake@params[["height_mm"]],
    units = "mm"
)

log_info("Writing between-haplotypes contrast table")
lapply(
    combinations(unique(age.metadata$Haplotype)),
    function(pair) {
      ages.1 <- age.metadata %>%
        filter(Haplotype == pair[1]) %>%
        drop_na(!!agecol) %>%
        pull(!!agecol)
      ages.2 <- age.metadata %>%
        filter(Haplotype == pair[2]) %>%
        drop_na(!!agecol) %>%
        pull(!!agecol)
      if (length(ages.1) == 0 || length(ages.2) == 0) {
        data.frame(
          `Haplotype 1` = pair[1], `Haplotype 2` = pair[2],
          Method = NA,
          Alternative = NA,
          Statistic = NA,
          `P-value` = NA
        )
      } else {
        wilcox.ages <- wilcox.test(ages.1, ages.2)
        data.frame(
          `Haplotype 1` = pair[1], `Haplotype 2` = pair[2],
          Method = wilcox.ages$method,
          Alternative = wilcox.ages$alternative,
          Statistic = wilcox.ages$statistic,
          `P-value` = wilcox.ages$p.value
        )
      }
    }
  ) %>%
  bind_rows %>%
  write_csv(snakemake@output[["report_haplotype_contrast"]])
