#!/usr/bin/env R

library(glue)
library(tidyverse)
library(ape)
library(logger)
library(ggtree)
library(castor)
library(phylobase)

TREE <- snakemake@input[["tree"]]
OUT.FILE <- snakemake@output[["tree_p4"]]

log_info("Reading Newick tree")
tree <- read.tree(TREE)

log_info("Correcting UShER tree edge lengths")
tree.corrected <- tree
tree.corrected$edge.length[is.na(tree.corrected$edge.length)] <- 0
tree.corrected$node.label <- NULL

log_info("Converting to phylo4")
tree.p4 <- as(tree.corrected, "phylo4")

log_info("Saving phylo4 tree")
save(tree.p4, file = OUT.FILE)

log_info("Finished")
