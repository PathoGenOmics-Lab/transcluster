#!/usr/bin/env python3


import json

import pandas as pd
import numpy as np
from scipy import stats as st


if __name__ == "__main__":

    # Read input files
    try:
        clusters = pd.read_csv(snakemake.input.clusters)
    except pd.errors.EmptyDataError:
        clusters = pd.DataFrame(columns=["cluster_id", "label"])
    extracted_ids = set(pd.read_csv(snakemake.input.extracted_ids)["original_id"])
    with open(snakemake.input.ids_file) as f:
        input_ids = set(line.strip() for line in f)

    # Print ID counts
    print(f"{len(extracted_ids)} extracted IDs")
    print(f"{len(input_ids)} total IDs")

    # Compose summary
    summary = {}
    ## Count observations
    summary["total_observations"] = len(input_ids)
    summary["analyzed_observations"] = len(extracted_ids)
    transmitted_ids = set(clusters["label"])
    non_transmitted_ids = extracted_ids - transmitted_ids
    summary["transmitted_observations"] = len(transmitted_ids)
    ## Count transmission events
    summary["n_clusters"] = clusters["cluster_id"].unique().size
    summary["n_emergences"] = summary["n_clusters"] + len(non_transmitted_ids)
    ## Calculate events size stats
    cluster_sizes = clusters["cluster_id"].value_counts().values
    if clusters["cluster_id"].unique().size != 0:
        summary["mean_cluster_size"] = float(np.mean(cluster_sizes))
        summary["sd_cluster_size"] = float(np.std(cluster_sizes))
        summary["median_cluster_size"] = float(np.median(cluster_sizes))
        summary["mode_cluster_size"] = st.mode(cluster_sizes, keepdims=False).mode.item()
    else:
        summary["mean_cluster_size"] = None
        summary["sd_cluster_size"] = None
        summary["median_cluster_size"] = None
        summary["mode_cluster_size"] = None

    # Save summary
    with open(snakemake.output.table, "w") as fw:
        json.dump(summary, fw)
