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
    
    # load all IGX CLuster datasets and save in dictionary
    igx_dict = {}
    igx_dict["V_80"] = pd.read_csv(data_path + "/IGX/Cluster_240425_V_80.tsv", sep="\t",
         low_memory=False)
    igx_dict["V_J_70"] = pd.read_csv(data_path + "/IGX/Cluster_240425_V_J_70.tsv", sep="\t",
         low_memory=False)
    igx_dict["V_70"] = pd.read_csv(data_path + "/IGX/Cluster_240425_V_70.tsv", sep="\t",
         low_memory=False)
    
    # preprocess datasets
    for i in igx_dict.keys():
        igx_dict[i] = preprocess_IGX_data(igx_dict[i])
    
    # calculate random clustering rate for all cluster sizes
    for i in igx_dict.keys():
        
        igx = igx_dict[i]
        igx_cluster_sizes = igx["Unique Cluster Id"].value_counts().to_list()
        
        igx_random_clustering_rates = calculate_average_random_clustering_rate(igx_cluster_sizes, igx["Unique Clone Id"].to_list(), 
            cluster_comparison[["Clone_ID_A", "Clone_ID_B"]])
        igx_random_clustering_rates.to_csv(data_path + "/Random_clustering_rates/Random_clustering_rates_IGX_" + str(i) + ".csv", 
            index=False)
