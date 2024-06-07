#!/usr/bin/env R

library(glue)
library(tidyverse)
library(ape)
library(logger)
library(ggtree)
library(phylobase)

library(doParallel)


EXTRACTED.IDS <- snakemake@input[["extracted_ids"]]
TREE.P4  <- snakemake@input[["tree_p4"]]
OUTPUT.CLUSTERS  <- snakemake@output[["clusters"]]
OUTPUT.SUMMARY  <- snakemake@output[["summary"]]
MIN.PROP <- snakemake@params[["min_prop"]]
MIN.SIZE <- snakemake@params[["min_size"]]
LOG_EVERY_MINUTES <- snakemake@params[["log_every_minutes"]]
NCPU <- snakemake@threads

registerDoParallel(NCPU)

calculate.clade.stats <- function(tree.p4, node, targets) {
  # returns prop and size
  descendant.tips <- names(descendants(tree.p4, node, type = "tips"))
  n.tips <- length(descendant.tips)
  c(sum(descendant.tips %in% targets) / n.tips, n.tips)
}


ParallelCladeResult <- function(queuedNode = NULL, selectedNode = NULL) {
  # source: https://stackoverflow.com/a/44101674
  me <- list(
    queuedNodes = queuedNode,
    selectedNodes = selectedNode
  )
  class(me) <- append(class(me), "ParallelCladeResult")
  return(me)
}


create.log.every <- function(n) {
  last_run_time <- Sys.time()
  call_count <- 0
  total.n <- n
  function(minutes) {
    # Update call count
    call_count <<- call_count + 1
    # Calculate whether to run or not
    current_time <- Sys.time()
    time_diff <- difftime(current_time, last_run_time, units = "mins")
    if (time_diff >= minutes) {
      delay <- round(time_diff - minutes, 1)
      pct <- round(call_count / total.n, 1)
      log_info("[log delay = {delay} min] Analyzed {call_count} of {total.n} internal nodes ({pct} %)")
      # Update the last run time
      last_run_time <<- current_time
    }
  }
}


compute.clusters <- function(tree.p4, node, targets, min.size = 1, min.prop =  0.9) {

  log_info("Initializing queue and results")
  node.type <- nodeType(tree.p4)
  log.every <- create.log.every(length(which(node.type == "internal")))
  results <- c()
  .q <- c()

  # Check 'node' first
  log_info("Checking provided node {node}")
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

  log_info("Starting loop for the rest of nodes")
  while (length(.q) != 0) {

    # Dequeue current node
    node <- .q[1]
    .q <- .q[-1]

    # Calculate and iterate over children
    node.children <- children(tree.p4, node)
    selected.node.children <- node.children[node.type[node.children] != "tip"]
    log_debug("Discarding {length(node.children) - length(selected.node.children)} child tips of {node}")
    par.output <- lapply(
      node.children,
      function(child) {
        log_debug("Exploring {node}: calculating {child} stats")
        stats <- calculate.clade.stats(tree.p4, child, targets)  # prop and size
        par.result <- ParallelCladeResult()
        if ((stats[1] >= min.prop) && (stats[2] >= min.size)) {
          log_debug("{child} qualifies (prop=stats[1]}, size={stats[2]})")
          log_debug("{child} to results")
          par.result$selectedNode <- child
        } else {
          log_debug("{child} does not qualify (prop=stats[1]}, size={stats[2]})")
          if (stats[1] * stats[2] < min.size) {
            # There are not even enough target nodes in the subclade to be selected
            log_debug("{child} does not contain enough targets (prop={stats[1]} * size={stats[2]} < min.size={min.size})")
          } else {
            # Here some subclade could be selected
            log_debug("{child} contains {stats[1]*stats[2]} targets, enqueuing")
            par.result$queuedNode <- child
          }
        }
        par.result
      }
    )
    log.every(LOG_EVERY_MINUTES)
    results <- c(results, unlist(sapply(par.output, function(x) x$selectedNode)))
    .q <- c(.q, unlist(sapply(par.output, function(x) x$queuedNode)))
  }
  log_info("Finished cluster computing")
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
cluster.table %>% write_csv(OUTPUT.CLUSTERS)
cluster.table %>%
  count(cluster_id) %>%
  write_csv(OUTPUT.SUMMARY)

log_info("All done")
