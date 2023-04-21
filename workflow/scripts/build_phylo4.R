#!/usr/bin/env R

library(glue)
library(phylobase)
library(tidyverse)
library(ape)
library(logger)

TREE <- snakemake@input[["tree"]]
EXTRACTED.IDS <- snakemake@input[["extracted_ids"]]
OUT.FILE <- snakemake@output[["tree_p4"]]

log_info("Reading target IDs")
ids <- read_csv(EXTRACTED.IDS)$modified_id
log_info("Read {length(ids)} IDs")

log_info("Reading Newick tree")
tree <- read.tree(TREE)
log_info("Tree has {length(tree$tip.label)} tips")

log_info("{sum(ids %in% tree$tip.label)}/{length(ids)} IDs are tips in tree")

log_info("Correcting UShER tree edge lengths")
tree$edge.length[is.na(tree$edge.length)] <- 0
tree$node.label <- NULL

log_info("Calculating target MRCA")
ids.mrca <- getMRCA(tree, ids)
log_info("Trimming tree")
tree <- extract.clade(tree, ids.mrca, root.edge = 1, collapse.singles = FALSE)
log_info("Trimmed tree has {length(tree$tip.label)} tips")

log_info("Converting to phylo4")
tree.p4 <- as(tree, "phylo4")

log_info("Saving phylo4 tree")
save(tree.p4, file = OUT.FILE)

log_info("Finished")
