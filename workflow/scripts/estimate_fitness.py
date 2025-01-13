#!/usr/bin/env python
#SBATCH --job-name corr
#SBATCH --mem 32G
#SBATCH --qos short
#SBATCH -c 1

"""
For each haplotype, for each cluster: obtain a score as the quotient
of the size between the number of sequences deposited in gisaid with collection
date between the first and last case of the cluster in the country or countries
(adding) where the cluster is observed
"""

from typing import List, Callable
import logging

import pandas as pd


# Other settings
DAYS_PADDING = pd.Timedelta(days=int(snakemake.wildcards.days_padding))

# Cluster detection analysis columns (hard-coded in scripts)
ANALYSIS_HAPLOTYPE_COL  = "Haplotype"
ANALYSIS_CLUSTER_ID_COL = "cluster_id"


class CountryBound:
    def __init__(self, country: str, start: pd.Timestamp, end: pd.Timestamp):
        self.country = country
        self.start = start
        self.end = end
    def __repr__(self):
        return f"CountryBound({self.country}[{self.start.date()}, {self.end.date()}])"


def extract_country(string):
    return string.split(" / ")[1] if " / " in string else string


def calculate_background_bounds(cluster: pd.DataFrame) -> List[CountryBound]:
    bounds = []
    for country, subdf in cluster.groupby("Country"):
        bounds.append(CountryBound(country, subdf["Date"].min(), subdf["Date"].max()))
    return bounds


def count_background_samples_slicing(metadata: pd.DataFrame, bounds: List[CountryBound]) -> int:
    unique_identifiers = set()
    for bound in bounds:
        min_bg_date = bound.start - DAYS_PADDING
        max_bg_date = bound.end + DAYS_PADDING
        new_ids = metadata[
            (min_bg_date < metadata["Date"]) & \
            (metadata["Date"] < max_bg_date) & \
            (metadata["Country"] == bound.country)
        ][snakemake.params.metadata_id_column].unique()
        unique_identifiers.update(new_ids)
        logging.debug(f"{len(new_ids)} samples in {bound.country} from {min_bg_date.date()} to {max_bg_date.date()}")
        if len(new_ids) == 0:
            matches_dates = sum((min_bg_date < metadata["Date"]) & (metadata["Date"] < max_bg_date))
            matches_loc = sum(metadata["Country"] == bound.country)
            logging.warning(f"{matches_loc} samples in {bound.country} and {matches_dates} samples from {min_bg_date.date()} to {max_bg_date.date()}")
    count = len(unique_identifiers)
    logging.debug(f"{count} samples in {bounds} +/- {DAYS_PADDING.days}")
    return count


def count_background_samples_adding(metadata: pd.DataFrame, bounds: List[CountryBound]) -> int:
    min_bg_date = min(bound.start for bound in bounds) - DAYS_PADDING
    max_bg_date = max(bound.end   for bound in bounds) + DAYS_PADDING
    locations = [bound.country for bound in bounds]
    count = len(
        metadata[
            (min_bg_date < metadata["Date"]) & \
            (metadata["Date"] < max_bg_date) & \
            metadata["Country"].isin(locations)
        ]
    )
    logging.debug(f"{count} samples in {locations} from {min_bg_date.date()} to {max_bg_date.date()}")
    return count


def calculate_corrected_fitness(cluster: pd.DataFrame, metadata: pd.DataFrame, func: Callable) -> float:
    cluster_size = len(cluster)
    bounds = calculate_background_bounds(cluster)
    background_count = func(metadata, bounds)
    if background_count > 0:
        return cluster_size / background_count
    else:
        logging.warning(f"Cannot calculate fitness <{func.__name__}> in {cluster['Country'].dropna().unique()} from {cluster['Date'].min().date()} to {cluster['Date'].max().date()} (no background samples)")
        return pd.NA


def analyze_cluster_results(cluster_data: pd.DataFrame, full_metadata: pd.DataFrame) -> pd.DataFrame:
    results = {"haplotype": [], "cluster_id": [], "fitness_adding": [], "fitness_slicing": []}
    for (haplotype, cluster_id), cluster in cluster_data.groupby([ANALYSIS_HAPLOTYPE_COL, ANALYSIS_CLUSTER_ID_COL]):
        min_date = cluster["Date"].min()
        max_date = cluster["Date"].max()
        locations = cluster["Country"].dropna().unique()
        if pd.notna(min_date) and pd.notna(max_date) and len(locations) > 0:
            logging.debug(f"Calculating fitness for haplotype {haplotype}/{cluster_id} (size {len(cluster)})")
            fitness_adding  = calculate_corrected_fitness(cluster, full_metadata, count_background_samples_adding)
            if len(locations) == 1:
                fitness_slicing = fitness_adding
            else:
                fitness_slicing = calculate_corrected_fitness(cluster, full_metadata, count_background_samples_slicing)
            logging.debug(f"Fitness adding = {fitness_adding:g} for haplotype {haplotype}/{cluster_id}")
            logging.debug(f"Fitness slicing = {fitness_slicing:g} for haplotype {haplotype}/{cluster_id}")
        else:
            logging.warning(f"Cannot calculate fitness in {locations} from {min_date.date()} to {max_date.date()} (not enough metadata)")
            fitness_adding = pd.NA
            fitness_slicing = pd.NA
        results["haplotype"].append(haplotype)
        results["cluster_id"].append(cluster_id)
        results["fitness_adding"].append(fitness_adding)
        results["fitness_slicing"].append(fitness_slicing)
    return pd.DataFrame.from_dict(results)


def read_cluster_analysis(path: str, sep: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=sep, dtype={ANALYSIS_CLUSTER_ID_COL: "Int64"})
    # Keep samples within transmission clusters (non-null cluster ID)
    df = df[df[ANALYSIS_CLUSTER_ID_COL].notna()]
    # Get countries
    df["Country"] = df[snakemake.params.metadata_location_column].apply(extract_country)
    # Set dates with year and month to the middle of the month and remove dates with only year
    date_items = df[snakemake.params.metadata_date_column].str.split("-").str.len()
    df.loc[date_items == 1, snakemake.params.metadata_date_column] = pd.NA
    df.loc[date_items == 2, snakemake.params.metadata_date_column] = df[snakemake.params.metadata_date_column] + "-15"
    df["Date"] = pd.to_datetime(df[snakemake.params.metadata_date_column], errors="raise")
    return df


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)-8s | %(message)s", level=logging.DEBUG)

    # Read metadata
    logging.info("Reading haplotype metadata with cluster IDs")
    clusters = read_cluster_analysis(snakemake.input.haplotype_metadata, sep=",")

    logging.info("Reading full metadata")
    metadata_sep = "," if snakemake.input.full_metadata.endswith(".csv") else "\t"
    metadata = pd.read_csv(
        snakemake.input.full_metadata,
        sep=metadata_sep, usecols=[snakemake.params.metadata_id_column, snakemake.params.metadata_date_column, snakemake.params.metadata_location_column]
    )
    metadata["Country"] = metadata[snakemake.params.metadata_location_column].apply(extract_country)

    # Set dates with year and month to the middle of the month and remove dates with only year
    date_items = metadata[snakemake.params.metadata_date_column].str.split("-").str.len()
    metadata.loc[date_items == 1, snakemake.params.metadata_date_column] = pd.NA
    metadata.loc[date_items == 2, snakemake.params.metadata_date_column] = metadata[snakemake.params.metadata_date_column] + "-15"
    metadata["Date"] = pd.to_datetime(metadata[snakemake.params.metadata_date_column], errors="raise")

    # Analyze results
    logging.info("Estimating fitness")
    fitness = analyze_cluster_results(clusters, metadata)

    # Write results
    logging.info("Writing results")
    fitness.to_csv(snakemake.output.fitness, index=False)
