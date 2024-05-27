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

    space2_dict = {}
    cutoffs = ["1", "1_5", "2", "2_5"]
    
    # load and preprocess all SPACE2 CLuster datasets and save in dictionary
    for cutoff in cutoffs:
        space2_dict[cutoff] = load_SPACE2_data("SPACE2_clustering_240425_cutoff_" + cutoff, data_path)
        space2_dict[cutoff] = preprocess_SPACE2_data(space2_dict[cutoff])

    # calculate random clustering rate for all cluster sizes
    for i in space2_dict.keys():
        
        space2 = space2_dict[i]

        space2_cluster_sizes = space2["cluster_by_rmsd"].value_counts().to_list()
        
        space2_random_clustering_rates = calculate_average_random_clustering_rate(space2_cluster_sizes, space2["Clone_ID"].to_list(), 
             cluster_comparison[["Clone_ID_A", "Clone_ID_B"]])
        space2_random_clustering_rates.to_csv(data_path + "/Random_clustering_rates/Random_clustering_rates_SPACE2_" + \
             str(i) + ".csv", index=False)