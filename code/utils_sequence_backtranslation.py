# -*- coding: utf-8 -*-
"""
Created on Thu May 23 2024

@author: Katharina Waury (k.waury@vu.nl)
"""
import itertools
import numpy as np
import pandas as pd
import random

from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

from utils import calculate_nucleotide_mutation_rate, concatenate_gene_segments, remove_partial_codons, translate_to_amino_acid


def add_matched_gene_information(epitope_df, reference_df, seq_col, verbose=True):
    
    matched_all = []
    
    for i, row in epitope_df.iterrows():
        
        if verbose:
            print(row["Epitope_ID"])

        best_genes, best_alignments = find_best_gene_alignment(row[seq_col], reference_df)
        assert len(best_genes) == len(best_alignments)  

        for i in range(len(best_alignments)):

            fraction_matching, mutation_num = examine_best_alignments(best_alignments[i], verbose=verbose)
            mismatch_list = find_mismatches(best_alignments[i][0][1], best_alignments[i][0][0])
            
            # add results to dataframe
            matched = [row["Epitope_ID"], row[seq_col], best_genes[i], best_alignments[i][0][1], fraction_matching, 
                mutation_num, mismatch_list]
            matched_all.append(matched)  

    epitope_df_matched = pd.DataFrame(matched_all, columns=["Epitope_ID", seq_col, "Gene_name", "Gene_seq", 
        "Fraction_match", "Mutation_num", "Mismatch_list"])

    return epitope_df_matched


def calculate_likelihood_of_current_mutation_order(mutation_pos_dict, mutation_order, seq, triplet_pos, triplet1, substitution_matrix):
    
    temporary_seq = seq
    cumulative_likelihood = 0
#     cumulative_likelihood = 1

    # per mutation order check cumulative substitution likelihood
    for i in mutation_order:
        mutation_pos = list(mutation_pos_dict.keys())[i]
        mutation_res = list(mutation_pos_dict.values())[i]
        
        # calculate and add likelihood of mutation
        likelihood = calculate_substitution_likelihood(triplet1, mutation_pos, mutation_res, substitution_matrix)
        
        cumulative_likelihood += likelihood
#         cumulative_likelihood *= likelihood
        
        print("\t\t\tLikelihood of current mutation (position: %i, residue: %s) is %.4f (Cumulative likelihood: %.4f)" 
              %(mutation_pos, mutation_res, likelihood, cumulative_likelihood))
        
        temporary_seq = update_nucleotide_sequence_with_nucleotide(temporary_seq, mutation_res, triplet_pos, 
            mutation_pos)
    
    return cumulative_likelihood




def calculate_mutation_number(triplet1, triplet2):
    
    mutations_required = sum(res1 != res2 for res1, res2 in zip(triplet1, triplet2))
    
    return mutations_required




def calculate_substitution_likelihood(triplet1, mutation_pos, mutation_res, substitution_matrix):
    
    # filter for position-specific substitution score
    substitution_matrix_pos_specific = filter_substitution_matrix(substitution_matrix, mutation_pos)
    likelihood = substitution_matrix_pos_specific.loc[triplet1, mutation_res]
            
    return likelihood




def concatenate_sequences(v_gene_seq, cdr3_seq, j_gene_seq):

    full_seq = v_gene_seq + cdr3_seq + j_gene_seq
    
    return full_seq


def create_final_nucleotide_sequences(row, cdr3_col, ref):
    
    row["Antibody_nucleotide_seq_V"], row[cdr3_col] = move_cysteine_to_cdr(row["Antibody_nucleotide_seq_V"], row[cdr3_col])
    
    row["Full_nucleotide_seq"] = concatenate_sequences(row["Antibody_nucleotide_seq_V"], row[cdr3_col], 
        row["Antibody_nucleotide_seq_J"])
    
#     row["Full_nucleotide_seq_with_leader"] = create_leader_seq(row["Full_nucleotide_seq"], row["Gene_name_V"], ref)
    
    return row


def create_leader_seq(full_seq, v_gene, ref):
    
    leader_seq = ref[ref["gene_name"] == v_gene]["L1_nt"].values[0] + ref[ref["gene_name"] == v_gene]["L2_nt"].values[0]
    leader_seq = remove_partial_codons(leader_seq)
    
    full_seq_with_leader = leader_seq + full_seq
    
    return full_seq_with_leader


def examine_best_alignments(alignment, verbose=True):
    
    # divide score (i.e. number of matching positions) by lenght of antibody sequence
    fraction_matching = alignment[0][2]/len(alignment[0][1])*100
    mutation_num = len(alignment[0][1]) - alignment[0][2]
    
    if verbose:
        print(format_alignment(*alignment[0], full_sequences=True))
        print("Fraction of matching amino acids: %.2f%%" %fraction_matching)
        print("Number of mutated positions: %i" %mutation_num)
        print("---------------------")
    
    return fraction_matching, mutation_num



def filter_substitution_matrix(substitution_matrix, mutation_pos):
    
    if mutation_pos == 0:
        substitution_matrix_pos_specific = substitution_matrix[substitution_matrix["Position"] == "Start"]
        
    elif mutation_pos == 1:
        substitution_matrix_pos_specific = substitution_matrix[substitution_matrix["Position"] == "Center"]
        
    elif mutation_pos == 2:
        substitution_matrix_pos_specific = substitution_matrix[substitution_matrix["Position"] == "End"]
            
    return substitution_matrix_pos_specific


def find_cdr3(df, ab_seq_col, cdr3_col):
    
    df["CDR3_start"] = df[ab_seq_col].find(df[cdr3_col])
    # check if any CDR3 regions could not be matched
    if df["CDR3_start"] == -1:
        print(df["Epitope_ID"])
    df["CDR3_end"] = df["CDR3_start"] + len(df[cdr3_col])
    
    return df


def find_best_gene_alignment(ab_seq_fragment, ref_genes):
    
    best_genes = []
    best_alignments = [] 
    best_score = 0
        
    for i, row in ref_genes.iterrows():
        current_gene = row["gene_name"]
        # retrieve only best alignment for current pair, returns Alignment object
        current_alignment = pairwise2.align.localxs(ab_seq_fragment, row["Amino_acid_seq"], -2, -2, one_alignment_only=True)
        current_score = current_alignment[0][2]
        
        # if current score is higher than previous best score, overwrite best score and best alignment
        if best_score < current_score:
            best_genes = [current_gene]
            best_alignments = [current_alignment]
            best_score = current_score
            
        # if current score is as good as best score, add alignemnt to list of best alignments    
        elif best_score == current_score:
            best_genes.append(current_gene)
            best_alignments.append(current_alignment)
        
        # skip scores below best score
        else:
            continue
    
    return best_genes, best_alignments



def find_mismatches(seqA, seqB):
    
    align = pairwise2.align.localxs(seqA, seqB, -2, -2, one_alignment_only=True)
    seqA_aligned = align[0][0]
    seqB_aligned = align[0][1]

    mismatch_dict = {}

    for i in range(len(seqB_aligned)):
        if seqA_aligned[i] != seqB_aligned[i]:
            mismatch_dict[i] = (seqA_aligned[i], seqB_aligned[i])
    
    return mismatch_dict


def find_most_likely_mutation_order(triplet1, triplet2, triplet_pos, mutation_pos_dict, seq, substitution_matrix):
    
    # create all possible orders of mutations
    combinations = list(itertools.permutations(range(len(mutation_pos_dict)), len(mutation_pos_dict)))
    
    mutation_order_likelihood_list = []

    for mutation_order in combinations:
        
        current_likelihood = calculate_likelihood_of_current_mutation_order(mutation_pos_dict, mutation_order, seq, 
            triplet_pos, triplet1, substitution_matrix)
        
        print("\t\tCurrent mutation order %s with likelihood %.4f" %(mutation_order, current_likelihood))
        mutation_order_likelihood_list.append(current_likelihood)    
        
    max_mutation_order_likelihood = max(mutation_order_likelihood_list)
    print("\t\tHighest triplet likelihood: %.4f" %max_mutation_order_likelihood)
    
    return max_mutation_order_likelihood


def find_mutation_position(triplet1, triplet2):
    
    mutation_pos_dict = {}

    for i in range(len(triplet1)):
        if triplet1[i] != triplet2[i]:
            mutation_pos_dict[i] = triplet2[i]
            
    if len(mutation_pos_dict) == 1:
        mutation_pos = list(mutation_pos_dict.keys())[0]
        mutation_res = list(mutation_pos_dict.values())[0]
        return mutation_pos, mutation_res
    
    elif len(mutation_pos_dict) > 1:
        return mutation_pos_dict

    
def find_triplets_with_least_number_of_mutations(triplet1, amino_acid, amino_acid_to_triplets):
    
    mutation_dict = {}
    
    for triplet2 in amino_acid_to_triplets[amino_acid]:
        mutations_required = calculate_mutation_number(triplet1, triplet2)
        mutation_dict[triplet2] = mutations_required
        
    # find lowest required mutation number
    lowest_mutation_number = min(mutation_dict.values())

    # only keep triplets that require minimum number of mutations  
    minimum_mutation_dict = {key: value for key, value in mutation_dict.items() if value == lowest_mutation_number}
    
    minimum_mutation_list = list(minimum_mutation_dict.keys())
    
    return minimum_mutation_list, lowest_mutation_number


def get_all_mismatch_positions(df):

    mismatch_list_all = []

    for i, row in df.iterrows():
        
        for key in row["Mismatch_list"]:

            value = row["Mismatch_list"][key] 
            if value[1] == "-":
                continue
            else:
                mismatch_list_all.append(key)
    
    return mismatch_list_all


def get_V_J_gene_fragment(df, ab_seq_col):
    
    df["Antibody_seq_V"] = df[ab_seq_col][:df["CDR3_start"]]
    df["Antibody_seq_J"] = df[ab_seq_col][df["CDR3_end"]:]
    
    return df



def insert_triplet(seq, triplet2, pos):
    
    updated_seq = seq[:pos*3] + triplet2 + seq[pos*3:]
    
    return updated_seq


def keep_most_likely_sequence_reconstruction(df):
    
    df.sort_values(["Epitope_ID", "Likelihood"], ascending=False, inplace=True)
    df.drop_duplicates(subset=["Epitope_ID"], keep="first", inplace=True)
    
    return df



def move_cysteine_to_cdr(v_gene, cdr3):
    
    v_gene_updated = v_gene[:-3]
    cdr3_updated = v_gene[-3:] + cdr3
    
    assert len(v_gene_updated) == len(v_gene) - 3 
    assert len(cdr3_updated) == len(cdr3) + 3
    
    return v_gene_updated, cdr3_updated


def read_cdr3_backtranslation_file(file, cdr3_nseq_col):
    
    headers = []
    sequences = []

    with open(file, "r") as f:
        header = ""
        sequence = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    headers.append(int(header))
                    sequences.append(sequence)
                header = line[1:]
                sequence = ""
            else:
                sequence += line
        if header:
            headers.append(int(header))
            sequences.append(sequence)

    cdr3_backtranslated = pd.DataFrame(data={"Epitope_ID":headers, cdr3_nseq_col:sequences})
    
    return(cdr3_backtranslated)

def reconstruct_from_reference_gene(mismatches, ref_seq, amino_acid_to_triplets, substitution_matrix):
    
    seq = ref_seq
    
    print("Gene mismatches:", mismatches)
    
    if len(mismatches) == 0:
        
        return seq, 1, 0  
    
    else:

        likelihood_list = []
        nucleotide_mutation_number = 0
        
        for i in mismatches.keys():
            
            reference_amino_acid = mismatches[i][0]
            antibody_amino_acid = mismatches[i][1]
            triplet1 = ref_seq[i*3:i*3+3]
            
            if reference_amino_acid == "-":
                print("Gap in reference sequence")
                # insert nucleotide triplet into sequence
                triplet2 = amino_acid_to_triplets[antibody_amino_acid][0]
                print("Triplet to be inserted into sequence:", triplet2)
                seq = insert_triplet(seq, triplet2, i)
                continue
                
            if antibody_amino_acid == "-":
                print("Gap in antibody sequence. Removing corresponding triplet from sequence")
                # replace corresponding triplet in sequence with placeholder
                seq = update_nucleotide_sequence_with_triplet(seq, "XXX", i)
                continue
                
            minimum_mutation_list, lowest_mutation_number = find_triplets_with_least_number_of_mutations(triplet1, 
                antibody_amino_acid, amino_acid_to_triplets)
            print("Mutation at position %i: %s" %(i, mismatches[i]))
            print("Triplet in reference:", triplet1)
            print("Potential triplets in antibody sequence:", minimum_mutation_list, 
                  "with %i mutation" %lowest_mutation_number)
            
            # if only one triplet has the minimum number of mutations, we use this triplet
            if len(minimum_mutation_list) == 1:
                
                triplet2 = minimum_mutation_list[0]
                
                if lowest_mutation_number == 1:
                    
                    mutation_pos, mutation_res = find_mutation_position(triplet1, triplet2)
                    likelihood = calculate_substitution_likelihood(triplet1, mutation_pos, mutation_res, substitution_matrix)
                    
                elif lowest_mutation_number > 1:
                    
                    print("\tComparing mutation order substitution likelihood...")
                    mutation_pos_dict = find_mutation_position(triplet1, triplet2)
                    
                    likelihood = find_most_likely_mutation_order(triplet1, triplet2, i, mutation_pos_dict, 
                        seq, substitution_matrix)
                    
                print("\tMutation likelihood %.4f" %likelihood)    
                seq = update_nucleotide_sequence_with_triplet(seq, triplet2, i)
               
            # if more than one triplet is possible, we need to check the substitution table 
            elif len(minimum_mutation_list) > 1: 
                
                if lowest_mutation_number == 1:              
                
                    print("\tComparing nucleotide substitution likelihood...")
                    triplet2, likelihood = select_most_likely_triplet(triplet1, i, minimum_mutation_list, seq, 
                        substitution_matrix)
                
                elif lowest_mutation_number > 1:
                
                    print("\tComparing mutation order substitution likelihood...")
                    triplet2, likelihood = select_triplet_with_most_likely_mutation_order(triplet1, i, 
                        minimum_mutation_list, seq, substitution_matrix)

                print("\tMost likely triplet is %s with likelihood %.4f" %(triplet2, likelihood))  

                seq = update_nucleotide_sequence_with_triplet(seq, triplet2, i)
                       
            likelihood_list.append(likelihood)
            nucleotide_mutation_number += lowest_mutation_number
            
        # remove gap placeholders from sequence 
        seq = seq.replace("X", "")
        
        # if only gaps at the end are part of the mismatches, the likelihood list will be empty
        if len(likelihood_list) == 0:
            return seq, 1, 0
        
        print("Cumulative likelihood: %.4f" %np.mean(likelihood_list))
        print("Nucleotide mutation number: %i" %nucleotide_mutation_number)
        print("-------------------")        
        
        return seq, np.mean(likelihood_list), nucleotide_mutation_number      
    

def reconstruct_nucleotide_sequence(row, ref, col_name, amino_acid_to_triplets, substitution_matrix):
    
    print("Epitope:", row["Epitope_ID"])
    
    # set reference gene as nucleotide sequence for antibody
    nucleotide_seq = retrieve_reference_nucleotide_sequence(ref, row["Gene_name"])

    row[col_name], row["Likelihood"], row["Mutation_num_nucleotide"] = reconstruct_from_reference_gene(row["Mismatch_list"], 
        nucleotide_seq, amino_acid_to_triplets, substitution_matrix)

    return row

    



def retrieve_reference_nucleotide_sequence(ref_df, gene):
    
    nucleotide_seq = ref_df[ref_df["gene_name"] == gene]["Shortened_seq"].values[0]
    
    return nucleotide_seq


def select_most_likely_triplet(triplet1, triplet_pos, minimum_mutation_list, seq, substitution_matrix):
    
    substitution_likelihood_dict = {}
    
    for triplet2 in minimum_mutation_list:
        
#         print("\t\tComparing %s and %s" %(triplet1, triplet2))
        
        # find position of mutation
        mutation_pos, mutation_res = find_mutation_position(triplet1, triplet2)
        likelihood = calculate_substitution_likelihood(triplet1, mutation_pos, mutation_res, substitution_matrix)
        print("\t\tTriplet %s with substitution likelihood %.4f" %(triplet2, likelihood))
        
        substitution_likelihood_dict[likelihood] = triplet2
        
    # compare substitution likelihoods and return most likely triplet
    most_likely_triplet = substitution_likelihood_dict[max(substitution_likelihood_dict.keys())] 
    
    return most_likely_triplet, max(substitution_likelihood_dict.keys())


def select_triplet_with_most_likely_mutation_order(triplet1, triplet_pos, minimum_mutation_list, seq, substitution_matrix):
    
    substitution_likelihood_dict = {}
    
    for triplet2 in minimum_mutation_list:
        
        print("\t\tComparing %s and %s" %(triplet1, triplet2))
        
        mutation_pos_dict = find_mutation_position(triplet1, triplet2)
        max_mutation_order_likelihood = find_most_likely_mutation_order(triplet1, triplet2, triplet_pos, mutation_pos_dict, 
            seq, substitution_matrix)
    
        substitution_likelihood_dict[max_mutation_order_likelihood] = triplet2
            
    # compare substitution likelihoods
    most_likely_triplet = substitution_likelihood_dict[max(substitution_likelihood_dict.keys())]         
    
    return most_likely_triplet, max(substitution_likelihood_dict.keys())




def update_nucleotide_sequence_with_nucleotide(seq, mutation_res, triplet_pos, mutation_pos):

    mutation_pos_in_seq = triplet_pos*3 + mutation_pos  
    updated_seq = seq[:mutation_pos_in_seq] + mutation_res + seq[mutation_pos_in_seq+1:]

    return updated_seq


def update_nucleotide_sequence_with_triplet(seq, triplet2, pos):
    
    updated_seq = seq[:pos*3] + triplet2 + seq[pos*3+3:]
    
    return updated_seq