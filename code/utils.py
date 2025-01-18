# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 14:10:40 2023

@author: Katharina Waury (k.waury@vu.nl)
"""
import itertools
import Levenshtein as ls
import numpy as np
import pandas as pd
import random
import re

from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


amino_acid_to_triplets = {
    "F": ["TTT", "TTC"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "I": ["ATT", "ATC", "ATA"],
    "M": ["ATG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "Y": ["TAT", "TAC"],
    "H": ["CAT", "CAC"],
    "Q": ["CAA", "CAG"],
    "N": ["AAT", "AAC"],
    "K": ["AAA", "AAG"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "C": ["TGT", "TGC"],
    "W": ["TGG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "*": ["TAA", "TAG", "TGA"]
}


amino_acid_palette = {
    "A": "Hydrophobic", 
    "C": "Hydrophobic", 
    "D": "Negative",
    "E": "Negative", 
    "F": "Aromatic", 
    "G": "Hydrophobic", 
    "H": "Positive", 
    "I": "Hydrophobic", 
    "K": "Positive", 
    "L": "Hydrophobic", 
    "M": "Hydrophobic",
    "N": "Polar",
    "P": "Hydrophobic",
    "Q": "Polar",
    "R": "Positive",
    "S": "Polar", 
    "T": "Polar",
    "V": "Hydrophobic", 
    "W": "Aromatic",
    "Y": "Aromatic"
}


def add_clone_id(epitope_id, mapping):
    
    try:
        return mapping[mapping["Epitope_ID"] == epitope_id]["Clone_ID"].values[0]
    except IndexError:
        return None
    

def add_cluster(clone_id, df, clone_col, cluster_col):
    
    try:
        cluster_id = df[df[clone_col] == clone_id][cluster_col].values[0]
    except:
        print(clone_id)
    return cluster_id


def add_gene_match(df, cluster_df):
        
    clone_A = df["Clone_ID_A"]
    clone_B = df["Clone_ID_B"]
    
    v_gene_A = cluster_df[cluster_df["Unique Clone Id"] == clone_A]["V Gene"].values[0]
    v_gene_B = cluster_df[cluster_df["Unique Clone Id"] == clone_B]["V Gene"].values[0]
    j_gene_A = cluster_df[cluster_df["Unique Clone Id"] == clone_A]["J Gene"].values[0]
    j_gene_B = cluster_df[cluster_df["Unique Clone Id"] == clone_B]["J Gene"].values[0]
    
    df["V_gene_match"] = (v_gene_A == v_gene_B)
    df["J_gene_match"] = (j_gene_A == j_gene_B)
        
    return df


def calculate_average_epitope_size(df):

    return (df["Length_A"] + df["Length_B"]) / 2


def calculate_average_random_clustering_rate(cluster_sizes, antibodies, ab_pairs, n=1000):
    
    assert sum(cluster_sizes) == len(antibodies), "The sum of the cluster sizes must match the number of antibodies"
    
    clustering_rate_list = []
    clustered_count_list = []
    
    for i in range(n):
        
        random_clusters = perform_random_clustering(cluster_sizes, antibodies)
        
        clustered_count, clustering_rate = calculate_random_clustering_rate(random_clusters, ab_pairs)
    
        clustered_count_list.append(clustered_count)
        clustering_rate_list.append(clustering_rate)

    clustering_rate_df = pd.DataFrame(data={"Clustered_count":clustered_count_list, "Clustering_rate":clustering_rate_list})
    
    return clustering_rate_df


def calculate_lower_bound_estimate_of_corrected_SHM(row):
    
    # calculate receptor length without the CDR3 region
    row["Full Sequence Nucleotides Length corrected"] = row["Full Sequence Nucleotides Length"] - \
        row["CDR3 Nucleotides Length"]
    
    # calculate corrected SHM rate when not considering CDR3 region
    row["SHM corrected (lower estimate)"] = row["V & J Gene SHM"] * row["Full Sequence Nucleotides Length"] / \
        row["Full Sequence Nucleotides Length corrected"]
        
    return row


def calculate_nucleotide_mutation_rate(row, mutatio_num_cols, seq_cols):
        
    mutation_counter = 0
    
    for i in mutatio_num_cols:
        mutation_counter += row[i]
        
    seq_length = 0
    
    for i in seq_cols:
        seq_length += len(row[i])
    
    mutation_rate = mutation_counter / seq_length * 100
    
    return mutation_rate


def calculate_random_clustering_rate(random_clusters, ab_pairs):
    
    clustered_count = 0
    
    for i, row in ab_pairs.iterrows():
        

        if check_if_antibodies_are_in_same_cluster(row["Clone_ID_A"], row["Clone_ID_B"], random_clusters):
            clustered_count += 1
    
    clustering_rate = clustered_count / len(ab_pairs)
    
    return clustered_count, clustering_rate


def calculate_relative_mismatch_position(row):
    
    relative_mismatch_pos = []
    
    for i in row["Mismatch_position"]:
        relative_pos = i/row["CDRH3_length"]
        relative_mismatch_pos.append(relative_pos)
    
    return relative_mismatch_pos


def calculate_sequence_distances(row, ab_seqs):
    
    # many epitopes/antibodies are lost because only antibody pairs with a PDB in SAbDab for each antibody are contained 
    if check_if_epitopes_have_antibody_sequence(row["Epitope_ID_A"], row["Epitope_ID_B"], ab_seqs):
        cdrh3_A = get_cdr3(row["Epitope_ID_A"], ab_seqs)
        cdrh3_B = get_cdr3(row["Epitope_ID_B"], ab_seqs)
        ab_seq_A = get_heavy_chain_sequence(row["Epitope_ID_A"], ab_seqs)
        ab_seq_B = get_heavy_chain_sequence(row["Epitope_ID_B"], ab_seqs)

        # calculate (normalized) distance between CDRH3 regions
        distance = ls.distance(cdrh3_A, cdrh3_B)
        distance_norm = ls.ratio(cdrh3_A, cdrh3_B)
        distance_norm_ab = ls.ratio(ab_seq_A, ab_seq_B)
        results_list = [row["Epitope_ID_A"], row["Epitope_ID_B"], cdrh3_A, cdrh3_B, distance, distance_norm, 
            distance_norm_ab]
            
    else:
        return None
        
    return results_list


def calculate_upper_bound_estimate_of_corrected_SHM(row):
    
    # calculate sum of all CDR region lengths
    row["CDR Nucleotides Length"] = row[["CDR1 Nucleotides Length", "CDR2 Nucleotides Length", \
        "CDR3 Nucleotides Length"]].sum()
    
    # calculate CDR region length without the CDR3 region
    row["CDR Nucleotides Length corrected"] = row["CDR Nucleotides Length"] - row["CDR3 Nucleotides Length"]
       
    # calculate corrected SHM rate when not considering CDR3 region
    row["SHM corrected (upper estimate)"] = row["V & J Gene SHM"] * row["CDR Nucleotides Length"] / \
        row["CDR Nucleotides Length corrected"]
    
    return row

        
def check_for_wrong_pairs_in_clusters(antibody_list, epitope_clone_mapping, cluster_dict): 
        
    # generate all combinations of antibodies in the cluster
    antibody_combinations = list(itertools.combinations(antibody_list, 2))
    false_positive_pairs = [] 
    
    # check if both antibodies are part of the antibody pair dataset
    for comb in antibody_combinations:
        if (epitope_clone_mapping["Clone_ID"].str.contains(comb[0]).any()) & \
            (epitope_clone_mapping["Clone_ID"].str.contains(comb[1]).any()):
            
            epitope_id_A = epitope_clone_mapping[epitope_clone_mapping["Clone_ID"] == comb[0]]["Epitope_ID"].values[0]
            epitope_id_B = epitope_clone_mapping[epitope_clone_mapping["Clone_ID"] == comb[1]]["Epitope_ID"].values[0]
            
            cluster_A = get_cluster(epitope_id_A, cluster_dict)
            cluster_B = get_cluster(epitope_id_B, cluster_dict)
            
            if cluster_A != cluster_B:
                false_positive_pairs.append([epitope_id_A, epitope_id_B])
                
    return false_positive_pairs


def check_if_antibodies_are_in_same_cluster(antibody_A, antibody_B, random_clusters):
    
    clustered = any((antibody_A in cluster) and (antibody_B in cluster) for cluster in random_clusters)
    
    return clustered


def check_if_epitopes_have_antibody_sequence(epitope_A, epitope_B, epitope_pdb_df):
    
    return (epitope_A in epitope_pdb_df["Epitope_ID"].values) & (epitope_B in epitope_pdb_df["Epitope_ID"].values)


def concatenate_gene_segments(df):
    
    V_regions = ["FR1_nt", "CDR1_nt", "FR2_nt", "CDR2_nt", "FR3_nt", "CDR3_V_nt"]
    full_seq = ""
    
    if df["gene_type"] == "V":
        for region in V_regions:
            full_seq += df[region]
        df["Full_seq"] = full_seq
        
    elif df["gene_type"] == "J":
        # last residue (3 nulceotides) of CDRH3 region belong to J gene according to IMGT definitions
        full_seq += df["CDR3_J_nt"][-3:] + df["FR4_nt"] 
        df["Full_seq"] = full_seq
            
    else:
        print("Gene type neither V or J. Gene will be skipped.")
        df["Full_seq"] = None
    
    return df


def create_antibody_pair_dataframe(epitope_pairs_df, ab_seqs):
    
    results_list = epitope_pairs_df.apply(calculate_sequence_distances, args=[ab_seqs], axis=1)
    # drop rows that returned None
    results_list.dropna(inplace=True)
    
    antibody_pair_df = pd.DataFrame(results_list.tolist(), columns=["Epitope_ID_A", "Epitope_ID_B", "CDRH3_A", 
    "CDRH3_B", "Distance", "Distance_normalized", "Distance_normalized_AB_seq"])
    
    epitope_set = set(antibody_pair_df["Epitope_ID_A"]).union(set(antibody_pair_df["Epitope_ID_B"]))
    print("Number of antibody pairs before filtering:", len(antibody_pair_df))
    print("Number of unique antibodies before filtering:", len(epitope_set))
    
    antibody_pair_df_filtered = filter_antibody_pair_dataframe(antibody_pair_df)

    epitope_set_filtered = set(antibody_pair_df_filtered["Epitope_ID_A"]).union(set(antibody_pair_df_filtered["Epitope_ID_B"]))
    print("Number of antibody pairs:", len(antibody_pair_df_filtered))
    print("Number of unique antibodies:", len(epitope_set_filtered))
    
    return antibody_pair_df_filtered, epitope_set_filtered


def create_immune_builder_dataframe(df):
    
    counter = 0
    
    antibody_list = []

    for i in df["Unique Clone Id"].unique():

        current_antibody = df[df["Unique Clone Id"] == i]

        # check that both heavy and light chain are present
        chains = df[df["Unique Clone Id"] == i]["Chain"].unique()
        if (sorted(chains) == ["Heavy", "Kappa"]) | (sorted(chains) == ["Heavy", "Lambda"]) == False:
            print("Could not find both heavy and light chain for the antibody. Skipping %s" %i)
            continue

        else:
            heavy_chain = current_antibody[current_antibody["Chain"] == "Heavy"]["Receptor Amino Acids"].values[0]
            light_chain = current_antibody[current_antibody["Chain"].isin(["Kappa", "Lambda"])]["Receptor Amino Acids"].values[0]

            current_antibody = [i, heavy_chain, light_chain]
            antibody_list.append(current_antibody)
        
        immune_builder_df = pd.DataFrame(antibody_list, columns = ["Clone_ID", "H", "L"]) 

        counter += 1

    return immune_builder_df


def create_list_of_epitope_IDs(epitope_ids, epitope_comparison_df):
        
    # get epitope IDs and add to list
    epitope_ids.extend(epitope_comparison_df["Epitope_ID_A"])
    epitope_ids.extend(epitope_comparison_df["Epitope_ID_B"])
    
    # remove duplicate epitope IDs
    epitope_ids = np.unique(np.array(epitope_ids)).tolist()
    
    return epitope_ids


def create_list_of_series_entries(series):
        
    exploded_series = series.explode().dropna()
    pos_list = exploded_series.tolist()
    
    return pos_list


def create_pie_labels(counts, n=3):

    labels = [x for x in counts.index[0:n]]
    labels.extend([""]*(len(counts) - n))
        
    return labels


def create_sequence_dataframe(ig_seq_dict, chain):

    cdr3_dict = {}
    ab_seq_dict = {}
    
    # create dictionary of full chain sequence and CDR3 region
    for i in ig_seq_dict.keys():
        ab_seq, cdr3_seq = get_sequence(ig_seq_dict[i])
        cdr3_dict[i] = cdr3_seq
        ab_seq_dict[i] = ab_seq
        
    # create dataframes from dictionary
    if chain == "H":
        cdr3_df = pd.DataFrame.from_dict(cdr3_dict, orient="index").reset_index().rename(columns={"index":"PDB", 0:"CDRH3"})
        ab_seq_df = pd.DataFrame.from_dict(ab_seq_dict, orient="index").reset_index().rename(columns={"index":"PDB", 
            0:"Heavy_chain"})
        ab_seq_df = ab_seq_df.merge(cdr3_df, on="PDB")
        ab_seq_df["CDRH3"] = ab_seq_df["CDRH3"].apply(remove_gaps)
        ab_seq_df["Heavy_chain"] = ab_seq_df["Heavy_chain"].apply(remove_gaps)
        
    if chain == "L":
        cdr3_df = pd.DataFrame.from_dict(cdr3_dict, orient="index").reset_index().rename(columns={"index":"PDB", 0:"CDRL3"})
        ab_seq_df = pd.DataFrame.from_dict(ab_seq_dict, orient="index").reset_index().rename(columns={"index":"PDB", 
            0:"Light_chain"})
        ab_seq_df = ab_seq_df.merge(cdr3_df, on="PDB")
        ab_seq_df["CDRL3"] = ab_seq_df["CDRL3"].apply(remove_gaps)
        ab_seq_df["Light_chain"] = ab_seq_df["Light_chain"].apply(remove_gaps)
        
    return ab_seq_df


def extract_pdb_ids_from_fasta(file_path):
    pdb_ids = []
    with open(file_path, "r") as fasta_file:
        for line in fasta_file:
            # check for header lines
            if line.startswith(">"):
                match = re.search(r"(\w{4})", line)
                if match:
                    pdb_id = match.group(1)
                    pdb_ids.append(pdb_id)
    return pdb_ids


def filter_antibody_pair_dataframe(antibody_df):
        
    # drop antibody pairs that have entirely the same sequence
    antibody_df_filtered = antibody_df[antibody_df["Distance_normalized_AB_seq"] < 1]
    
    # remove identical CDRH3 pairs (might be differences in antibody sequence but still redundant)
    antibody_df_filtered = antibody_df_filtered.drop_duplicates(subset=["CDRH3_A", "CDRH3_B"], keep="first")
    
    return antibody_df_filtered 


def filter_epitopes_for_high_overlap(epitope_comparison_df, cutoff=0.9, how="Jaccard_score"):
    
    epitope_pairs_similar = epitope_comparison_df[epitope_comparison_df[how] > cutoff] 
    
    return epitope_pairs_similar


def filter_fastq(input_file, output_file, numbers):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        header = ""
        sequence = ""
        quality = ""

        for line in infile:
            if line.startswith("@"):
                header = line.strip()
                sequence = infile.readline().strip()
                plus_line = infile.readline().strip()
                quality = infile.readline().strip()
                        

def filter_for_antibodies_with_two_entries(sabdab_iedb):
    
    pdb_counts = sabdab_iedb["pdb"].value_counts()
    sabdab_iedb_two_entries = sabdab_iedb[sabdab_iedb["pdb"].isin(pdb_counts[pdb_counts == 2].index)]

    return sabdab_iedb_two_entries


def filter_for_same_class_antibodies(sabdab_iedb):

    # remove PDBs where the antibody information on chains does not match --> cannot be the same antibody
    keep = []
    print("Dropping:", end = " ")

    for pdb in sabdab_iedb["pdb"].unique():
        current_df = sabdab_iedb[sabdab_iedb["pdb"] == pdb]
        if (len(current_df["heavy_subclass"].unique()) == 1) & \
        (len(current_df["light_subclass"].unique()) == 1) & \
        (len(current_df["light_ctype"].unique()) == 1):
                keep.append(pdb)
        else:
            print(pdb, end = " ")

    # keep only the PDBs likely to contain same antibodies in two different entries
    sabdab_iedb_same_class = sabdab_iedb[sabdab_iedb["pdb"].isin(keep)]
    
    return sabdab_iedb_same_class


def filter_sequence_information(df_merged):
    
    # drop antibodies with same sequence
    df_merged.drop_duplicates(subset=["Heavy_chain", "Light_chain"], keep="first", inplace=True)
    
    # if epitope ID was mapped to multiple PDBs, keep only one
    df_merged = df_merged.drop_duplicates(subset=["Epitope_ID"], keep="first")
    
    # drop antibodies with non-standard amino acid
    df_merged = df_merged[~df_merged["Heavy_chain"].str.contains("X")]
    df_merged = df_merged[~df_merged["Light_chain"].str.contains("X")]
        
    return df_merged


def filter_SAbDab_entries(sabdab):
    
    # only paired H & L chains and defined antigen chains
    sabdab_filtered = sabdab.dropna(subset=["Hchain", "Lchain", "antigen_type"])
    # only protein antigens
    sabdab_filtered = sabdab_filtered[sabdab_filtered["antigen_type"] == "protein"]
    # no single-chain variable fragments
    sabdab_filtered = sabdab_filtered[sabdab_filtered["scfv"] == False]
    # only X-ray crystal structures
    sabdab_filtered = sabdab_filtered[sabdab_filtered["method"] == "X-RAY DIFFRACTION"]
    
    # # minimum resolution
    # sabdab_filtered = sabdab_filtered[sabdab_filtered["resolution"] != "NOT"]
    # sabdab_filtered["resolution"] = sabdab_filtered["resolution"].astype("float64")
    # sabdab_filtered = sabdab_filtered[sabdab_filtered["resolution"] < 3.5]

    return sabdab_filtered


def filter_SAbDab_entries_by_number(sabdab):
    
    # count occurrences of each pdb value
    pdb_counts = sabdab["pdb"].value_counts()
    # filter rows with one or two occurrences of pdb
    sabdab_filtered_by_numbers = sabdab[sabdab["pdb"].isin(pdb_counts[pdb_counts <= 2].index)]
    # one row indicates that only one chain assignment exists
    # two rows indicate that there are differing chain assignments between authors and PDB 
    # or that two antibodies are present in structure   
    
    # only one or two occurences of PDBs left
    assert all(sabdab_filtered_by_numbers["pdb"].value_counts() <= 2)
    
    return sabdab_filtered_by_numbers


def find_clone_ID_of_substring(main_string, df):
    
    substring_list = df["Receptor Amino Acids"]
    matching_seq = [substring in main_string for substring in substring_list]
    
    clone_id = df[matching_seq]["Unique Clone Id"].values[0]
    
    return clone_id   


def find_PDB(epitope_id, epitope_pdb_mapping):
    
    return epitope_pdb_mapping[epitope_pdb_mapping["Epitope_ID"] == epitope_id]["PDB"].values[0]

    
def find_sequence_cutoff(mismatch_list):
    
    # create list of mismatch position sorted in reverse
    mismatch_positions = list(mismatch_list)
    mismatch_positions.sort(reverse=True)
    
    # update cutoff position based on number of gaps at the sequence end
    for i in mismatch_positions:
        aligned_res = mismatch_list[i][1]
        if aligned_res == "-":
            cutoff = i
            
        else:
            break
    
    return cutoff


def get_cdr3(epitope, epitope_df):
    
    return epitope_df[epitope_df["Epitope_ID"] == epitope]["CDRH3"].values[0]


def get_cluster(epitope, cluster_dict):
    
    for key, value in cluster_dict.items():
        if epitope in value:
            return key


def get_heavy_chain_sequence(epitope, epitope_df):
    
    return epitope_df[epitope_df["Epitope_ID"] == epitope]["Heavy_chain"].values[0]


def get_IEDB_ID(url):
    
    idx = url.rsplit("/", maxsplit=1)[1]
    
    return int(idx)


def get_PDB(epitope, antibody_seq_df):
    
    return antibody_seq_df[antibody_seq_df["Epitope_ID"] == epitope]["PDB"].values[0]


def get_sequence(pos_res_list):
    
    # IMGT numbering
    cdr3_pos = ["105", "106",  "107",  "108",  "109",  "110",  "111", "111A", "111B", "111C", "111D", "111E", "111F", "112F",
                "112E", "112D", "112C", "112B",  "112A",  "112", "113", "114", "115", "116", "117"]
    
    cdr3_seq = ""
    ab_seq = ""
    
    for i in pos_res_list:
        ab_seq += i[1]
        if i[0] in cdr3_pos:
            cdr3_seq += i[1]
    
    return ab_seq, cdr3_seq 


def is_any_substring_in_list(main_string, substring_list):
    
    return any(substring in main_string for substring in substring_list)


def load_SPACE2_data(filename, data_path):
    
    df = pd.read_csv(data_path + filename + ".csv")
    df.rename(columns={"ID":"Clone_ID"}, inplace=True)
    
    return df


def merge_sequence_information(epitope_pdb, heavy_chain_df, light_chain_df):

    df_merged = epitope_pdb.merge(heavy_chain_df, on="PDB", how="inner")
    df_merged = df_merged.merge(light_chain_df, on="PDB", how="inner")
    print("Number of ANARCI numbered antibody sequences:", len(df_merged))
    
    df_filtered = filter_sequence_information(df_merged)
    
    print("Number of filtered ANARCI numbered antibody sequences:", len(df_filtered))
    
    return df_filtered


def perform_random_clustering(cluster_sizes, antibodies):
    
    random_clusters = []
    remaining_antibodies = list(antibodies)

    for size in cluster_sizes:
        
        random_cluster = random.sample(remaining_antibodies, size)
        random_clusters.append(random_cluster)
        remaining_antibodies = [antibody for antibody in remaining_antibodies if antibody not in random_cluster]
        
    assert len(random_clusters) == len(cluster_sizes)

    return random_clusters


def preprocess_IGX_data(df):
    
    germline, cluster = separate_germline_seq(df) 
    cluster_heavy = cluster[cluster["Chain"] == "Heavy"]
    
    return cluster_heavy


def preprocess_SPACE2_data(df):
    
    # create integer mapping for unique strings
    integer_mapping = {string: i for i, string in enumerate(df["cluster_by_rmsd"].unique())}
    # replace strings with integers 
    df["cluster_by_rmsd"] = df["cluster_by_rmsd"].map(integer_mapping)

    # extract ID from path string
    df["Clone_ID"] = df["Clone_ID"].apply(rename_ID)
    
    return df


def print_p_val(p_val):
 
    if p_val < 0.0001:
        return "< 0.0001"
    else:
        return "%.4f" % p_val


def read_anarci_results_file(file_path, chain):

    ig_seq_dict = {}
    pos_list = []
    res_list = []

    with open(file_path, "r") as f:
        for line in f.readlines():
            # save zipped list of position and residue to dictionary
            if "#" in line:
                if len(pos_list) != 0:
                    combined_list = list(zip(pos_list, res_list))
                    ig_seq_dict[current_pdb] = combined_list
                
                # reset lists
                pos_list = []
                res_list = []

                # get PDB ID from first line
                if line[6:7] == "_":
                    current_pdb = line[2:6]
                # skip other information lines
                else:
                    continue

            elif line[:1] == chain:
                line_condensed = line.replace(" ", "").strip()
                pos_list.append(line_condensed[1:-1])
                res_list.append(line_condensed[-1:]) 
                
    # Add the last entry to the dictionary if the file doesn't end with "#"
    if len(pos_list) != 0:
        combined_list = list(zip(pos_list, res_list))
        ig_seq_dict[current_pdb] = combined_list
    
    return ig_seq_dict


def read_reference_file(file_path):

    genes = []
    sequences = []
    
    with open(file_path, "r") as fasta_file:
        
        gene = ""  
        sequence = "" 
        
        for line in fasta_file:
            line = line.strip() 
            if line.startswith(">"):
                if gene:
                    genes.append(gene)
                    sequences.append(sequence)
                gene = line[1:]  
                sequence = ""
            else:
                sequence += line
        
        # add last entry to lists
        if gene:
            genes.append(gene)
            sequences.append(sequence)

    # create dataframe from the dictionary
    ref = pd.DataFrame({"gene_name": genes, "Full_seq": sequences})
    
    return ref


def remove_gaps(string):
    
    return string.replace("-", "")   


def remove_partial_codons(string):

    # check for complete codons (1 codon = 3 nucleotides)
    if len(string) % 3 == 0:
        return string
    else:
        # remove partial codon
        shortened_string = string[:-(len(string) % 3)]
        assert len(shortened_string) % 3 == 0
        return shortened_string
    
    
def remove_sequence_end(seq, cutoff):
    
    new_seq = seq[:cutoff]

    return new_seq


def remove_sequence_redundancy(df, seq_col):
    
    cutoff = find_sequence_cutoff(df["Mismatch_list"])

    df["Gene_seq"] = remove_sequence_end(df["Gene_seq"], cutoff)
    df[seq_col] = remove_sequence_end(df[seq_col], cutoff*3)
    
    return df


def rename_ID(string):
    
    # remove file path
    _, new_string = string.rsplit("/", 1)
    # remove pdb extension
    new_string, _ = new_string.split(".", 1)
    
    return new_string


def retrieve_clone_ID(epitope_id, df):
    
    clone_id = df[df["Epitope_ID"] == epitope_id]["Clone_ID"].values[0]
    
    return clone_id 


def retrieve_protein_seq(epitope_id, df):
    
    seq = df[df["Epitope_ID"] == epitope_id]["Protein_Seq"].values[0]
    
    return seq 


def retrieve_substitution_matrix_score(fivemer, mutation_res, substitution_matrix):
        
    likelihood = substitution_matrix.loc[fivemer, mutation_res]
    
    return likelihood


def sequence_translation(df):
    
    df = df.apply(concatenate_gene_segments, axis=1)
    df["Shortened_seq"] = df["Full_seq"].apply(remove_partial_codons)
    df["Amino_acid_seq"] = df["Shortened_seq"].apply(translate_to_amino_acid)
    
    # check if beginning and end of amino acid sequence are as expected
    df["Amino_acid_seq_start"] = df["Amino_acid_seq"].str[:3]
    df["Amino_acid_seq_end"] = df["Amino_acid_seq"].str[-3:]
    
    return df


def separate_germline_seq(df):
    germline = df[df["Germline Sequence"] == True]
    cluster = df[df["Germline Sequence"] == False]

    print("Germline sequences:", len(germline))
    print("Unique antibodies:", len(cluster["Unique Clone Id"].unique()))
    
    return germline, cluster


def translate_to_amino_acid(nt_seq):
    
    # create Seq object from Biopython module
    dna = Seq(nt_seq)
    # return sequence string, not Seq object
    aa_seq = str(dna.translate(stop_symbol="_"))
    
    return aa_seq