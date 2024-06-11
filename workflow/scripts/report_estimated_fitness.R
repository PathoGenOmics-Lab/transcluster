library(tidyverse)
library(logger)

Sys.setlocale("LC_TIME", "English")
theme_set(theme_classic())

combinations <- function(elements) {
  # Generate all possible pairs of elements
  pairs <- combn(elements, 2)
  combos <- list()
  for (i in seq_len(ncol(pairs))) {
    combos[[i]] <- pairs[, i]
  }
  return(combos)
}

log_info("Reading data")
estimated.fitness <- read_csv(snakemake@input[["fitness"]])

log_info("Plotting")
estimated.fitness %>%
  rename(Adding = fitness_adding, Slicing = fitness_slicing) %>%
  pivot_longer(cols = c(Adding, Slicing)) %>%
  ggplot(aes(Haplotype, value, color = name)) +
    geom_boxplot() +
    scale_y_log10() +
    annotation_logticks(sides = "l") +
    scale_x_discrete(drop = FALSE) +
    ylab("Estimated transmission fitness") +
    # 1) K-W
    stat_compare_means(method = "wilcox", paired = FALSE) +
    # 2) Pairwise
    stat_compare_means(
      comparisons = combinations(unique(estimated.fitness$Haplotype)),
      method = "wilcox", paired = FALSE
    )

log_info("Writing report")
ggsave(
  snakemake@output[["report"]],
  width = snakemake@params[["width_mm"]],
  height = snakemake@params[["height_mm"]],
  units = "mm"
)

log_info("Writing summary")
estimated.fitness %>%
  group_by(Haplotype) %>%
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
