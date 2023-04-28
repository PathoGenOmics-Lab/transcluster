#!/usr/bin/env python3


import json
from pathlib import Path
import pandas as pd
import numpy as np


if __name__ == "__main__":

    # Read input files
    clusters_path = Path(snakemake.input.cluster_dir) / "clusters.csv"
    clusters = pd.read_csv(clusters_path)
    extracted_ids = set(pd.read_csv(snakemake.input.extracted_ids)["original_id"])
    with open(snakemake.input.ids_file) as f:
        input_ids = set(line.strip() for line in f)

    # Print ID counts
    print(f"{len(extracted_ids)} extracted IDs")
    print(f"{len(input_ids)} extracted IDs")

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
    summary["mean_cluster_size"] = float(np.mean(cluster_sizes))
    summary["sd_cluster_size"] = float(np.std(cluster_sizes))
    summary["median_cluster_size"] = float(np.median(cluster_sizes))
    summary["mode_cluster_size"] = int(np.max(cluster_sizes))

    # Save summary
    with open(snakemake.output.table, "w") as fw:
        json.dump(summary, fw)
