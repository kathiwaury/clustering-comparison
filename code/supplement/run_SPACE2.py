# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 16:09:40 2024

@author: Katharina Waury (k.waury@vu.nl)
"""
# adapted from https://github.com/oxpig/SPACE2

import glob
import os
import pandas as pd
import SPACE2
import time


if __name__ == "__main__":

    timer_list = []
    
    start = time.time()

    antibody_models = glob.glob(os.path.dirname(os.getcwd()) + "/data/ImmuneBuilder_R2/*.pdb")
    results_df = SPACE2.agglomerative_clustering(antibody_models, cutoff=2)
    
    results_df.to_csv(os.path.dirname(os.getcwd()) + "/data/SPACE2/SPACE2_clustering_R2_cutoff_2.csv", index=False)
    
    end = time.time()
    time_passed = end - start
    timer_list.append(time_passed)
    print("All structures clustered. SPACE2 took %.3f seconds to run" % time_passed, flush=True)





