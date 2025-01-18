# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 15:11:55 2024

@author: Katharina Waury (k.waury@vu.nl)
"""
import os
import pandas as pd
import sys

sys.path.append("..")
from utils import *


if __name__ == "__main__":
    
    data_path = os.path.dirname(os.path.dirname(os.getcwd())) + "/data"

    # load all data
    # cluster comparison
    cluster_comparison = pd.read_csv(data_path + "/antibody_pairs/Antibody_pairs_cluster_comparison.csv")
    
    # IGX Cluster
    cluster_full = pd.read_csv(data_path + "/clustering/IGX_Cluster_clustering.tsv", sep="\t", low_memory=False)
    germline, cluster = separate_germline_seq(cluster_full)
    cluster_heavy = cluster[cluster["Chain"] == "Heavy"]
    
    # SAAB+
    saabplus = pd.read_csv(data_path + "/clustering/SAABplus_clustering_heavy_chains.tsv", sep="\t", index_col=0)
    saabplus_clustered = saabplus[~saabplus["Clusters"].isna()] # in BAZIS "None" is interpreted as None, not a string
    
    # SPACE2
    space2 = pd.read_csv(data_path + "/clustering/SPACE2_clustering.csv")
    space2.rename(columns={"ID":"Clone_ID"}, inplace=True)

    # preprocess data sets
    clones_seq = cluster_heavy[["Unique Clone Id", "Receptor Amino Acids"]]
    clones_seq.columns = ["Clone_ID", "Protein_Seq"]
    saabplus_clustered_clones = clones_seq.merge(saabplus_clustered, on="Protein_Seq")
    
    integer_mapping = {string: i for i, string in enumerate(space2["cluster_by_rmsd"].unique())}
    space2["cluster_by_rmsd"] = space2["cluster_by_rmsd"].map(integer_mapping)
    space2["Clone_ID"] = space2["Clone_ID"].apply(rename_ID)

    # create list of cluster sizes
    igx_cluster_sizes = cluster_heavy["Unique Cluster Id"].value_counts().to_list()
    saabplus_cluster_sizes = saabplus_clustered_clones["Clusters"].value_counts().to_list()
    space2_cluster_sizes = space2["cluster_by_rmsd"].value_counts().to_list()
    
    # calculate random clustering rates
    # IGX Cluster
    igx_random_clustering_rates = calculate_average_random_clustering_rate(igx_cluster_sizes, 
        cluster_heavy["Unique Clone Id"].to_list(), cluster_comparison[["Clone_ID_A", "Clone_ID_B"]])
    igx_random_clustering_rates.to_csv(data_path + "/random_clustering_rates/Random_clustering_rates_IGX.csv", index=False)

    # SAAB+
    saabplus_random_clustering_rates = calculate_average_random_clustering_rate(saabplus_cluster_sizes, 
    saabplus_clustered_clones["Clone_ID"].to_list(), cluster_comparison[["Clone_ID_A", "Clone_ID_B"]])
    saabplus_random_clustering_rates.to_csv(data_path + "/random_clustering_rates/Random_clustering_rates_SAABplus.csv", index=False)

    # SPACE2
    space2_random_clustering_rates = calculate_average_random_clustering_rate(space2_cluster_sizes, 
        space2["Clone_ID"].to_list(), cluster_comparison[["Clone_ID_A", "Clone_ID_B"]])
    space2_random_clustering_rates.to_csv(data_path + "/random_clustering_rates/Random_clustering_rates_SPACE2.csv", index=False)

