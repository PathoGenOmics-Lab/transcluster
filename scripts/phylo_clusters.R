#!/usr/bin/env R

library(glue)
library(tidyverse)
library(ape)
library(logger)
library(ggtree)
library(phylobase)


IDS.FILE <- snakemake@input[["ids_file"]]
TREE.P4  <- snakemake@input[["tree_p4"]]
OUT.DIR  <- snakemake@output[["out_dir"]]
MIN.PROP <- snakemake@params[["min_prop"]]
MIN.SIZE <- snakemake@params[["min_size"]]

LOG.EVERY.SECONDS <- 1


calculate.clade.stats <- function(tree.p4, node, targets) {
  # returns prop and size
  descendant.tips <- names(descendants(tree.p4, node, type = "tips"))
  n.tips <- length(descendant.tips)
  c(sum(descendant.tips %in% targets) / n.tips, n.tips)
}


compute.clusters <- function(tree.p4, node, targets, min.size = 1, min.prop =  0.9) {

  log_debug("Initializing queue, visited nodes and results")
  all.nodes <- nodeId(tree.p4, type = "all")
  .visited <- lapply(all.nodes,
                     function(x) FALSE)
  results <- c()
  print(.visited)
  n.nodes <- length(all.nodes)
  t <- Sys.time()

  # Check 'node' first
  log_debug("Checking provided node {node}")
  stats <- calculate.clade.stats(tree.p4, node, targets)
  if ((stats[1] >= min.prop) & (stats[2] >= min.size)) {
    log_debug("Subclades of node {node} do qualify (prop={round(stats[1], 2)}, size={stats[2]})")
    descs <- descendants(tree.p4, node, type = "all")
    log_debug("Marking {descs} as visited")
    .visited[descs] <- TRUE
    results <- c(results, node)
  } else {
    log_debug("Subclades of {node} do NOT qualify (prop={round(stats[1], 2)}, size={stats[2]})")
    log_debug("Marking {node} as visited")
    .visited[node] <- TRUE
  }

  log_debug("Starting loop for the rest of nodes")
  .q <- c(node)
  while (any(.q)) {
    # Log every LOG.EVERY.SECONDS seconds
    if (difftime(Sys.time(), t, units = "secs") > LOG.EVERY.SECONDS) {
      current.n.visited <- length(.visited)
      pct.done <- 100 * current.n.visited / n.nodes
      log_info("Visited {current.n.visited} nodes, search {round(pct.done, 2)}% complete")
      t <- Sys.time()
    }

    # Select current node and dequeue
    node <- .q[length(.q)]
    .q <- .q[-length(.q)]

    # Calculate and iterate over children
    node.children <- children(tree.p4, node)
    for (child in node.children) {
      if (!.visited[[child]]) {
        # Node is not visited
        stats <- calculate.clade.stats(tree.p4, child, targets)
        if ((stats[1] >= min.prop) & (stats[2] >= min.size)) {
          log_debug("Subclades of node {child} do qualify (prop={round(stats[1], 2)}, size={stats[2]})")
          # Mark as visited all members of every subclade with enough targets
          # (this way we can avoid exploring any subclade)
          descs <- descendants(tree.p4, child, type = "all")
          log_debug("Marking {descs} as visited")
          .visited[descs] <- TRUE
          # Add node to results
          results <- c(results, child)
        } else {
          log_debug("Subclades of {child} do NOT qualify (prop={round(stats[1], 2)}, size={stats[2]})")
          log_debug("Marking {child} as visited")
          # Enqueue node for further exploration and mark as visited
          .q <- c(.q, child)
          .visited[child] <- TRUE
        }
      } else {
        log_debug("Skipping node {child}")
      }
    }
  }
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


log_threshold(INFO)

log_info("Starting")

log_info("Creating output directory")
dir.create(OUT.DIR)

log_info("Reading Phylo4 tree")
load(TREE.P4)

log_info("Reading IDs")
targets <- read_lines(IDS.FILE)

log_info("Finding starting node")
target.mrca <- MRCA(tree.p4, targets)

log_info("Calculating clusters")
cluster.nodes <- compute.clusters(tree.p4, target.mrca, targets, MIN.SIZE, MIN.PROP)
cluster.table <- create.cluster.table(tree.p4, cluster.nodes)

log_info("Writing clusters")
cluster.table %>% write_csv(glue("{OUT.DIR}/clusters.csv"))
cluster.table %>%
  count(cluster_id) %>%
  write_csv(glue("{OUT.DIR}/cluster_summary.csv"))

log_info("All done")
