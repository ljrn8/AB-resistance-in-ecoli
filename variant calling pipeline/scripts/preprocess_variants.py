import numpy as np
import pandas as pd
import os
import pickle
import heapq
import pickle


def get_hashmap():
    snp_hash = {}

    for count, file in enumerate(os.scandir("../dataset")):
        genome_id, variant_type = file.name.rstrip(".txt").split("_")
        
        if variant_type != 'snps':
            continue
        
        with open("../dataset/" + file.name, "r") as f:
            
            for row in f:

                # TODO numpy call
                position = None 

                if position in snp_hash:
                    continue 
                else:
                    snp_hash[position] = np.array(1508)
            



def view_matrix(snp_matrix):
    i = 0
    for key in snp_matrix:
        variants = snp_matrix[key]
        print(f'{key}: {len(variants)}, {variants[:20]}')
        
        if i == 10:
            return
        i += 1


def generate_snp_matrix():
    with open("../results/positions_sorted", "rb") as loci_file:
        loci = pickle.load(loci_file)
    
    print(np.array(loci), len(loci))
    
    snp_matrix = {}
    tokenize = {
        b'A': 1,
        b'T': 2,
        b'C': 3,
        b'G': 4,
    }
    
    for position in loci:
        # unsigned byte
        snp_matrix[position] = np.zeros(1508, dtype=np.dtype('B'))
    
    print("finished snp matrix prep;")
    view_matrix(snp_matrix)
    
    for count, file in enumerate(os.scandir("../dataset")):
        genome_id, variant_type = file.name.rstrip(".txt").split("_")
        
        if variant_type != 'snps':
            continue
        
        if count % 1 == 0:
            print(f"\n\t[{count/2}/1508]")
            view_matrix(snp_matrix)
        
        with open("../dataset/" + file.name, "r") as f:
            X = np.loadtxt(f, dtype=int, usecols=(0, 2),
                converters={
                    2: lambda x: tokenize[x]
                }
            )
            # hash for snp
            for position, snp in X:
                snp_matrix[position][count] = tokenize[snp]
                
                
    with open("../results/snp_matrix", "wb") as snp_matrix_file:
        pickle.dump(snp_matrix, snp_matrix_file)
        
        
        
        # TODO issues for when you return
        # variant count problem
        # line 107 gets 2 nucleotides at once to tokinze, makes no sense
