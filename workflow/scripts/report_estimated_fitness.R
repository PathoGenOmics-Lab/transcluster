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

log_info("Reading data")
estimated.fitness <- lapply(
  snakemake@input[["fitnesses"]],
  function(path) {
    read_csv(path, col_types = c(
      "haplotype" = "c", "cluster_id" = "f",
      "fitness_adding" = "d", "fitness_slicing" = "d"
      )
    )
  }
) %>% bind_rows
n.haplotypes <- length(snakemake@input[["fitnesses"]])

log_info("Formatting plot data")
plot.data <- estimated.fitness %>%
  rename(Adding = fitness_adding, Slicing = fitness_slicing,
         Haplotype = haplotype) %>%
  pivot_longer(cols = c(Adding, Slicing))

log_info("Building global reports")
if (nrow(plot.data) != 0) {
  p.simple <- plot.data %>%
    ggplot(aes(Haplotype, value, color = name)) +
      geom_boxplot() +
      scale_y_log10() +
      annotation_logticks(sides = "l") +
      scale_x_discrete(drop = FALSE) +
      ylab("Estimated transmission fitness") +
      # 1) Compare between fitness types
      stat_compare_means(method = "wilcox", paired = FALSE)

  p.report <- p.simple +
      # 2) Compare pairwise
      stat_compare_means(
        comparisons = combinations(unique(estimated.fitness$haplotype)),
        method = "wilcox", paired = FALSE
      )

  p.slicing <- plot.data %>%
    filter(name == "Slicing") %>%
    ggplot(aes(Haplotype, value)) +
      geom_boxplot() +
      scale_y_log10() +
      annotation_logticks(sides = "l") +
      scale_x_discrete(drop = FALSE) +
      ylab("Estimated transmission fitness") +
      # Compare pairwise
      stat_compare_means(
        comparisons = combinations(unique(estimated.fitness$haplotype)),
        method = "wilcox", paired = FALSE
      )

  p.adding <- plot.data %>%
    filter(name == "Adding") %>%
    ggplot(aes(Haplotype, value)) +
      geom_boxplot() +
      scale_y_log10() +
      annotation_logticks(sides = "l") +
      scale_x_discrete(drop = FALSE) +
      ylab("Estimated transmission fitness") +
      # Compare pairwise
      stat_compare_means(
        comparisons = combinations(unique(estimated.fitness$haplotype)),
        method = "wilcox", paired = FALSE
      )

} else {
  log_warn("Plot data is empty")
  p.simple <- plot.data %>% ggplot(aes(Haplotype, value))
  p.report <- p.simple
  p.slicing <- p.simple
  p.adding <- p.simple
}

log_info("Writing global reports")
ggsave(
  snakemake@output[["report_simple"]],
  plot = p.simple,
  width = n.haplotypes * snakemake@params[["width_per_haplotype_mm"]],
  height = snakemake@params[["height_mm"]],
  units = "mm"
)

ggsave(
  snakemake@output[["report"]],
  plot = p.report,
  width = n.haplotypes * snakemake@params[["width_per_haplotype_mm"]],
  height = snakemake@params[["height_mm"]],
  units = "mm"
)

log_info("Writing fitness-adding report")
ggsave(
  snakemake@output[["report_adding"]],
  plot = p.adding,
  width = n.haplotypes * snakemake@params[["width_per_haplotype_mm"]],
  height = snakemake@params[["height_mm"]],
  units = "mm"
)

log_info("Writing fitness-slicing report")
ggsave(
  snakemake@output[["report_slicing"]],
  plot = p.slicing,
  width = n.haplotypes * snakemake@params[["width_per_haplotype_mm"]],
  height = snakemake@params[["height_mm"]],
  units = "mm"
)

log_info("Writing global report data")
plot.data %>% write_csv(snakemake@output[["report_data"]])

log_info("Writing summary")
estimated.fitness %>%
  group_by(haplotype) %>%
  summarize(
    # Fitness adding
    mean_fitness_adding = mean(fitness_adding, na.rm = TRUE),
    median_fitness_adding = median(fitness_adding, na.rm = TRUE),
    gmean_fitness_adding = exp(mean(log(fitness_adding), na.rm = TRUE)),
    sd_fitness_adding = sd(fitness_adding, na.rm = TRUE),
    sem_fitness_adding = sd_fitness_adding / sqrt(length(which(!is.na(fitness_adding)))),
    # Fitness slicing
    mean_fitness_slicing = mean(fitness_slicing, na.rm = TRUE),
    median_fitness_slicing = median(fitness_slicing, na.rm = TRUE),
    gmean_fitness_slicing = exp(mean(log(fitness_slicing), na.rm = TRUE)),
    sd_fitness_slicing = sd(fitness_slicing, na.rm = TRUE),
    sem_fitness_slicing = sd_fitness_slicing / sqrt(length(which(!is.na(fitness_slicing))))
  ) %>%
  write_csv(snakemake@output[["summary"]])

log_info("Writing contrast table")
lapply(
    combinations(unique(estimated.fitness$haplotype)),
    function(pair) {
      # Slicing
      slicing.1 <- estimated.fitness %>%
        filter(haplotype == pair[1]) %>%
        drop_na(fitness_slicing) %>%
        pull(fitness_slicing)
      slicing.2 <- estimated.fitness %>%
        filter(haplotype == pair[2]) %>%
        drop_na(fitness_slicing) %>%
        pull(fitness_slicing)
      if (length(slicing.1) == 0 || length(slicing.2) == 0) {
        wilcox.slicing <- list(method = NA, alternative = NA, statistic = NA, p.value = NA)
      } else {
        wilcox.slicing <- wilcox.test(slicing.1, slicing.2)
      }
      # Adding
      adding.1 <- estimated.fitness %>%
        filter(haplotype == pair[1]) %>%
        drop_na(fitness_adding) %>%
        pull(fitness_adding)
      adding.2 <- estimated.fitness %>%
        filter(haplotype == pair[2]) %>%
        drop_na(fitness_adding) %>%
        pull(fitness_adding)
      if (length(slicing.1) == 0 || length(slicing.2) == 0) {
        wilcox.adding <- list(method = NA, alternative = NA, statistic = NA, p.value = NA)
      } else {
        wilcox.adding <- wilcox.test(adding.1, adding.2)
      }
      # Final dataframe
      rbind(
          # Slicing
          data.frame(
            `Haplotype 1` = pair[1], `Haplotype 2` = pair[2],
            Type = "slicing",
            Method = wilcox.slicing$method,
            Alternative = wilcox.slicing$alternative,
            Statistic = wilcox.slicing$statistic,
            `P-value` = wilcox.slicing$p.value
          ),
          # Adding
          data.frame(
            `Haplotype 1` = pair[1], `Haplotype 2` = pair[2],
            Type = "adding",
            Method = wilcox.adding$method,
            Alternative = wilcox.adding$alternative,
            Statistic = wilcox.adding$statistic,
            `P-value` = wilcox.adding$p.value
          )
      )
    }
  ) %>%
  bind_rows %>%
  write_csv(snakemake@output[["report_contrast"]])
