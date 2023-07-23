import numpy as np
import pandas as pd
import os
import pickle
import heapq
import pickle

def get_loci(only_snps=True):
    
    # hashable check ? 
    unique_set = set()
    loci = []

    # do 1 pass first 
    with open("../dataset/0_snps.txt", "r") as f:
        loci = list(np.loadtxt(f, dtype=int, usecols=(0,)))
        print('should look like this: ', np.array(loci))
        
    for count, file in enumerate(os.scandir("../dataset")):
        genome_id, variant_type = file.name.rstrip(".txt").split("_")
        
        if only_snps and variant_type != 'snps':
            continue
        
        if count % 20 == 0:
            # only reflects snp count
            print(f'\n{count/2}/1508: ', loci[:10], loci[-10:-1], len(loci))
        
        with open("../dataset/" + file.name, "r") as f:
        
            burst = list(np.loadtxt(f, dtype=int, usecols=(0,)))

            # extract only unqiue loci
            common_set = unique_set.intersection(burst)
            unique_snps = list(set(burst) - common_set)                 
            
            print(len(unique_snps), ", ", end="")
                
            for snp in unique_snps:
                unique_set.add(snp)
                
            loci.extend(unique_snps)

            # Convert the extended array into a min-heap
            heapq.heapify(loci)
            
    with open("../results/positions", "wb") as loci_file:
        pickle.dump(loci, loci_file)
        


def sort_heap(heap_file="../results/loci"):
    with open(heap_file, "rb") as loci_file:
        loci_heap = pickle.load(loci_file)
        
    loci = [heapq.heappop(loci_heap) for i in range(len(loci_heap))]
    with open("../results/positions_sorted", "wb") as loci_file:
        pickle.dump(loci, loci_file)
    

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
            
            X = list(np.loadtxt(f, dtype=int, usecols=(0, 2),
                converters={
                    2: lambda x: tokenize[x]
                }
            ))

            # hash for snp
            for position, snp in X:
                snp_matrix[position][count] = tokenize[snp]
                
                
    with open("../results/snp_matrix", "wb") as snp_matrix_file:
        pickle.dump(snp_matrix, snp_matrix_file)
        
        
        
        # TODO issues for when you return
        # variant count problem
        # line 107 gets 2 nucleotides at once to tokinze, makes no sense