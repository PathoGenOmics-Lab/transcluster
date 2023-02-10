#!/usr/bin/env R

library(glue)
library(tidyverse)
library(ape)
library(logger)
library(ggtree)
library(phylobase)


IDS.DIR  <- snakemake@input[["ids_dir"]]
TREE.P4  <- snakemake@input[["tree_p4"]]
OUT.DIR  <- snakemake@output[["out_dir"]]
MIN.PROP <- snakemake@params[["min_prop"]]
MIN.SIZE <- snakemake@params[["min_size"]]
DATASET  <- snakemake@wildcards[["dataset"]]


calculate.clade.stats <- function(tree.p4, node, targets) {
  # returns prop and size
  descendant.tips <- names(descendants(tree.p4, node, type = "tips"))
  n.tips <- length(descendant.tips)
  c(sum(descendant.tips %in% targets) / n.tips, n.tips)
}


compute.clusters <- function(tree.p4, node, targets, min.size = 1, min.prop =  0.9) {

  log_debug("Initializing queue, visited nodes and results")
  .q <- c(node)
  .visited <- c(node, names(tipLabels(tree.p4)))  # add tips to visited
  results <- c()

  while (! identical(.q, numeric(0))) {
    # Select current node and dequeue
    node <- .q[length(.q)]
    .q <- .q[-length(.q)]

    # Calculate and iterate over children
    node.children <- children(tree.p4, node)
    for (child in node.children) {
      # If node is not visited and is internal
      if (!child %in% .visited) {
        stats <- calculate.clade.stats(tree.p4, child, targets)
        if ((stats[1] >= min.prop) & (stats[2] >= min.size)) {
          log_debug("Subclades of node {child} do qualify (prop={round(stats[1], 2)}, size={stats[2]})")
          # Mark as visited all members of every subclade with enough targets
          # (this way we can avoid exploring any subclade)
          log_debug("Marking {descendants(tree.p4, child, type = 'all')} as visited")
          .visited <- c(.visited, descendants(tree.p4, child, type = "all"))
          # Add node to results
          results <- c(results, child)
        } else {
          log_debug("Subclades of {child} do NOT qualify (prop={round(stats[1], 2)}, size={stats[2]})")
          # Enqueue node for further exploration and mark as visited
          .q <- c(.q, child)
          .visited <- c(.visited, child)
        }
      } else {
        log_debug("Skipping node {child}")
      }
    }
  }
  results
}


create.cluster.table <- function(tree.p4, cluster.nodes) {
  df <- data.frame()
  i <- 1
  for (node in cluster.nodes) {
    labels <- names(descendants(tree.p4, node, type = "tips"))
    df[labels, "label"] <- labels
    df[labels, "cluster_id"] <- i
    i <- i + 1
  }
  df
}


log_threshold(INFO)

log_info("Starting")

log_info("Creating output directory")
dir.create(OUT.DIR)

log_info("Reading Phylo4 tree")
load(TREE.P4)

log_info("Reading haplotype IDs")
haplotypes <- data.frame()
for (file.name in list.files(IDS.DIR, full.names = TRUE)) {
    # Read
    haplotype.name <- tools::file_path_sans_ext(basename(file.name))
    labels <- read_lines(file.name)
    # Build table for ggtree
    haplotypes[labels, "node"] <- labels
    haplotypes[labels, "Haplotype"] <- haplotype.name
}

log_info("Finding tree root")
tree.root <- rootNode(tree.p4)

for (haplotype.name in unique(haplotypes$Haplotype)) {
  log_info("Calculating clusters for {haplotype.name}")
  cluster.nodes <- compute.clusters(tree.p4, node, targets, 1, 0.9)
  log_info("Writing clusters for {haplotype.name}")
  create.cluster.table(tree.p4, cluster.nodes) %>% write_csv(glue("{OUT.DIR}/{haplotype.name}.csv"))
}

log_info("All done")
