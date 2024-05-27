# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:38:28 2023

@author: Katharina Waury (k.waury@vu.nl)
"""

#!/usr/bin/env python
# coding: utf-8

# import libraries
import numpy as np
import os
import pandas as pd
import re


def epitope_to_dict(epitope_desc):
    
    epitope_dict = {}
    
    # remove all spaces as these are not always consistent
    epitope_desc = epitope_desc.replace(" ", "")
    # split on commas to get single residues of epitope
    epitope_list = epitope_desc.split(sep=",")
    
    for res in epitope_list:
        
        # check if epitope was described in correct manner
        try:
            assert res[0] in "ACDEFGHIKLMNPQRSTUVWY"
            int(res[1:])
        except:
            print("The epitope %s was not described in a correct manner and was skipped" %(epitope_desc))
            return None
        
        # add residue amino acids and position to the epitope dict
        epitope_dict[int(res[1:])] = res[0]
    
    return epitope_dict
    

def filter_epitopes_for_antigen_subset(epitope_df, antigen, organism):
    
    epitope_df_filtered = epitope_df[(epitope_df["Epitope - Molecule Parent"] == antigen) & 
        (epitope_df["Epitope - Species"] == organism)]
    
    print("Number of epitopes:", len(epitope_df_filtered))
    
    return epitope_df_filtered


def calculate_shared_residues(epitope_dict_A, epitope_dict_B):
    
    shared_res = {k: epitope_dict_A[k] for k in epitope_dict_A if k in epitope_dict_B and 
                epitope_dict_A[k] == epitope_dict_B[k]}
    
    return shared_res


def calculate_jaccard_score(epitope_dict_A, epitope_dict_B):
    
    intersection = len(set(epitope_dict_A.items()).intersection(set(epitope_dict_B.items())))
    union = len(set(epitope_dict_A.items()).union(set(epitope_dict_B.items())))
    
    return intersection / union


def create_results_series(epitope_df, i, j):

    epitope_a = epitope_df.iloc[i]["Epitope - Dict"]
    epitope_b = epitope_df.iloc[j]["Epitope - Dict"]

    shared_res = calculate_shared_residues(epitope_a, epitope_b)
    jaccard_score = calculate_jaccard_score(epitope_a, epitope_b)

    results_list = [epitope_df.iloc[i]["Epitope ID - IEDB IRI"], epitope_df.iloc[j]["Epitope ID - IEDB IRI"], epitope_a, epitope_b,
                    len(epitope_a), len(epitope_b), shared_res, len(shared_res),
                    len(shared_res)/len(epitope_a), len(shared_res)/len(epitope_b),
                    np.mean([len(shared_res)/len(epitope_a), len(shared_res)/len(epitope_b)]),
                    jaccard_score]

    return results_list


if __name__ == "__main__":
    
    # import IEDB data
    epitopes = pd.read_csv(os.path.dirname(os.path.dirname(os.getcwd())) + "/data/IEDB/IEDB_epitopes_discontinuous_240215.csv")

    # set how many top antigens should be checked
    top = 50
    
    # identify highly annotated antigens
    top_antigens = epitopes[["Epitope - Molecule Parent", "Epitope - Species"]].value_counts()[:top]

    # define columns for result dataframe
    columns = ["Antigen", "Species", "Epitope_ID_A", "Epitope_ID_B", "Epitope_A", "Epitope_B", "Length_A", "Length_B", 
               "Shared_residues", "Overlap", "Fraction_overlap_A", "Fraction_overlap_B", "Fraction_overlap_avg", "Jaccard_score"]
    
    for antigen in top_antigens.index:
        
        print(antigen[0], antigen[1])
        epitope_df_filtered = filter_epitopes_for_antigen_subset(epitopes, antigen[0], antigen[1])[["Epitope ID - IEDB IRI", 
            "Epitope - Name"]]
        
        epitope_df_filtered["Epitope - Dict"] = epitope_df_filtered["Epitope - Name"].apply(epitope_to_dict)
        # remove incorrect epitopes
        epitope_df_filtered.dropna(subset=["Epitope - Dict"], inplace=True)

        epitope_comparison_list = [] # to save pairwise comparison results to

        # pairwise comparison of all epitope residues
        for i in range(len(epitope_df_filtered)):
            for j in range(len(epitope_df_filtered)):
                if j > i:
                    # get pairwise epitope similarity
                    results_series = create_results_series(epitope_df_filtered, i, j)
                    results_series = [antigen[0], antigen[1]] + results_series
                    epitope_comparison_list.append(results_series)

        # create results dataframe from list of results series
        epitope_comparison_df = pd.DataFrame(epitope_comparison_list, columns=columns)

        # save dataframes to CSV
        molecule = re.sub("[ ,/,(,)]", "", antigen[0])
        species = antigen[1].replace(" ", "")
        epitope_comparison_df.to_csv(os.path.dirname(os.path.dirname(os.getcwd())) + 
            "/data/IEDB/Epitope_comparison/Epitope_comparison_" + molecule + "_" + species + ".csv", index=False)
