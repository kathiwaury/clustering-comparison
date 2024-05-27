# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 10:10:26 2024

@author: Katharina Waury (k.waury@vu.nl)
"""
# adapted from https://github.com/oxpig/ImmuneBuilder

import numpy as np
import os
import pandas as pd
import time

from ImmuneBuilder import ABodyBuilder2


def create_output_filename(clone_id):
    
    output_file = clone_id + ".pdb"
    
    return output_file


def build_sequence_dict(heavy_chain, light_chain):
    
    sequences = {"H" : heavy_chain, "L" : light_chain}
    
    return sequences
    

def predict_antibody_structure(sequences):

    try:
    	antibody = predictor.predict(sequences)
    except:
        return None
    
    return antibody


if __name__ == "__main__":

    timer_list = []
    
    predictor = ABodyBuilder2()
    
    immune_builder_df = pd.read_csv(os.path.dirname(os.getcwd()) + "/data/SPACE2/ImmuneBuilder_input.csv")
    
    print("Number of clones:", len(immune_builder_df[4400:]), flush=True)
    
    for i, row in immune_builder_df[4400:].iterrows():

        if i % 1000 == 0:
            print("Number of predicted clones:", i, flush=True)
        
#        start = time.time()
        
        sequences = build_sequence_dict(row["H"], row["L"])
        antibody_structure = predict_antibody_structure(sequences)

        if antibody_structure is None:
            print("Failed structure prediction for %s" %row["Clone_ID"])
        else:
            antibody_structure.save(os.path.dirname(os.getcwd()) + "/data/ImmuneBuilder_R2/" + row["Clone_ID"] + ".pdb")
        
#        end = time.time()
#        time_passed = end - start
#        timer_list.append(time_passed)
#        print("Current clone %s took %.3f seconds to run" % (i, time_passed), flush=True)

# print("Average time required per structure prediction: %.3f seconds" % np.mean(timer_list))

print("All structures predicted.", flush=True)
