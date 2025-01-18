# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 17:09:52 2025

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
    
    # load all IGX cluster datasets and save in dictionary
    igx = pd.read_csv(data_path + "/clustering/IGX_Cluster_clustering_80.tsv", sep="\t", low_memory=False)
    
    # preprocess IGX dataset
    igx = preprocess_IGX_data(igx)
    igx.rename(columns={"Unique Clone Id":"Clone_ID"}, inplace=True)
    
    space2 = load_SPACE2_data("/clustering/SPACE2_clustering", data_path)
    space2 = preprocess_SPACE2_data(space2)
    
    igx_space2 = igx.merge(space2, on="Clone_ID")
    
    print(igx_space2[:5])    
    
    # calculate random clustering rate for only CDR length
    space2_cluster_sizes = space2["cluster_by_length"].value_counts().to_list()
    
    random_clustering_rates = calculate_average_random_clustering_rate(space2_cluster_sizes, space2["Clone_ID"].to_list(), 
         cluster_comparison[["Clone_ID_A", "Clone_ID_B"]])
    random_clustering_rates.to_csv(data_path + "/random_clustering_rates/Random_clustering_rates_SPACE2_length.csv", index=False)
    
    
    # create combined (CDRH3 sequence identity + CDR lengths) clustering column
    igx_space2["Cluster_ID_combined"] = igx_space2["Unique Cluster Id"].astype(str) + "_" + igx_space2["cluster_by_length"]
    print(len(igx_space2["Cluster_ID_combined"].unique()))
    
    clusters = igx_space2.groupby("Cluster_ID_combined")
    print(len(clusters))
    
    # calculate random clustering rate for all cluster sizes
    combined_cluster_sizes = igx_space2["Cluster_ID_combined"].value_counts().to_list()
        
    random_clustering_rates = calculate_average_random_clustering_rate(combined_cluster_sizes, igx_space2["Clone_ID"].to_list(), 
        cluster_comparison[["Clone_ID_A", "Clone_ID_B"]])
    random_clustering_rates.to_csv(data_path + "/random_clustering_rates/Random_clustering_rates_IGX_80_SPACE2_length.csv", 
        index=False)


