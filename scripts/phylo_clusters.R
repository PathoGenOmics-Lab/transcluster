#!/usr/bin/env R

library(glue)
library(tidyverse)
library(ape)
library(logger)
library(ggtree)
library(phylobase)


EXTRACTED.IDS <- snakemake@input[["extracted_ids"]]
TREE.P4  <- snakemake@input[["tree_p4"]]
OUT.DIR  <- snakemake@output[["out_dir"]]
MIN.PROP <- snakemake@params[["min_prop"]]
MIN.SIZE <- snakemake@params[["min_size"]]


calculate.clade.stats <- function(tree.p4, node, targets) {
  # returns prop and size
  descendant.tips <- names(descendants(tree.p4, node, type = "tips"))
  n.tips <- length(descendant.tips)
  c(sum(descendant.tips %in% targets) / n.tips, n.tips)
}


compute.clusters <- function(tree.p4, node, targets, min.size = 1, min.prop =  0.9) {

  log_debug("Initializing queue and results")
  node.type <- nodeType(tree.p4)
  results <- c()
  .q <- c()

  # Check 'node' first
  log_debug("Checking provided node {node}")
  stats <- calculate.clade.stats(tree.p4, node, targets)
  if ((stats[1] >= min.prop) && (stats[2] >= min.size)) {
    log_debug("Node {node} DOES qualify (prop={round(stats[1], 2)}, size={stats[2]})")
    log_debug("Adding {node} to results")
    results <- c(results, node)
  } else if (stats[1]*stats[2] < min.size) {
    log_debug("Subclades of node {node} do not contain enough targets")
  } else {
    log_debug("Node {node} does NOT qualify (prop={round(stats[1], 2)}, size={stats[2]})")
    log_debug("Enqueuing {node}")
    .q <- c(node)
  }

  log_debug("Starting loop for the rest of nodes")
  while (length(.q) != 0) {

    # Dequeue current node
    node <- .q[1]
    .q <- .q[-1]

    # Calculate and iterate over children
    node.children <- children(tree.p4, node)
    for (child in node.children) {

      if (node.type[child] == "tip") {
        log_debug("Exploring {node}: skipping tip {child}")
        next
      }

      log_debug("Exploring {node}: calculating {child} stats")
      stats <- calculate.clade.stats(tree.p4, child, targets)  # prop and size

      if ((stats[1] >= min.prop) && (stats[2] >= min.size)) {
        log_debug("{child} qualifies (prop={round(stats[1], 6)}, size={stats[2]})")
        log_debug("{child} to results")
        results <- c(results, child)
      } else {
        log_debug("{child} does not qualify (prop={round(stats[1], 6)}, size={stats[2]})")

        if (stats[1] * stats[2] < min.size) {
          # There are not even enough target nodes in the subclade to be selected
          log_debug("{child} does not contain enough targets ({stats[1]*stats[2]})")
        } else {
          # Here some subclade could be selected
          log_debug("{child} contains {stats[1]*stats[2]} targets")
          log_debug("Enqueuing {child}")
          .q <- c(.q, child)
        }
      }
    }
  }
  log_debug("Finished cluster computing")
  results
}


create.cluster.table <- function(tree.p4, cluster.nodes) {
  df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(df) <- c("label", "cluster_id")
  i <- 1
  for (node in cluster.nodes) {
    labels <- names(descendants(tree.p4, node, type = "tips"))
    df[labels, "label"] <- labels
    df[labels, "cluster_id"] <- i
    i <- i + 1
  }
  df
}


log_threshold(DEBUG)

log_info("Starting")

log_info("Creating output directory")
dir.create(OUT.DIR)

log_info("Reading Phylo4 tree")
load(TREE.P4)

log_info("Reading target IDs")
extracted.ids <- read_csv(EXTRACTED.IDS)
targets <- extracted.ids$modified_id
log_info("{length(targets)} targets read")

log_info("Finding root node")
tree.root <- rootNode(tree.p4)
log_info("Root node: {tree.root}")

log_info("Calculating clusters")
cluster.nodes <- compute.clusters(tree.p4, tree.root, targets, MIN.SIZE, MIN.PROP)

log_info("Generating table")
cluster.table <- create.cluster.table(tree.p4, cluster.nodes)

if (length(cluster.nodes) == 0) {
  log_info("No clusters found")
  cluster.table <- cbind(cluster.table, modified_id = character(0))
} else {
  log_info("{length(unique(cluster.table$cluster_id))} clusters found")
  cluster.table <- cluster.table %>%
    left_join(extracted.ids, by = c("label" = "modified_id"))
}

log_info("Writing clusters")
cluster.table %>% write_csv(glue("{OUT.DIR}/clusters.csv"))
cluster.table %>%
  count(cluster_id) %>%
  write_csv(glue("{OUT.DIR}/cluster_summary.csv"))

log_info("All done")
