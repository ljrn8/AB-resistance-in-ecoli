import numpy as np
import pandas as pd
import os
import pickle

# try one with snps and compare with snps + indels ? 

# try unsigned np.int32
# snp_matrix = np.array(range(0, 1509), dtype=np.int64)
# print(snp_matrix)

# df = pd.DataFrame(columns=['genome_id', 'snps']) 
# df['genome_id'] = pd.Series(range(0, 1509))
# print(df)

class loci_stack(list):
    def insertion_sort(self, item) -> None:
        self.append(item)
        for step in range(1, len(self)):
            key = self[step]
            j = step - 1
            while j >= 0 and key < self[j]:
                self[j + 1] = self[j]
                j = j - 1
            self[j + 1] = key 






def add_snps() -> np.array:
    
    snp_matrix = {}
    
    for count, file in enumerate(os.scandir("../dataset")):
        genome_id, variant_type = file.name.rstrip(".txt").split("_")
        
        if variant_type != 'snps':
            continue
        
        if count % 100 == 0:
            print(f'dataframe at file count {count} = ', df.info(), df)
        
        with open("../dataset/" + file.name, "r") as f:
        
            for line in f:
                position, reference, alternate, quality = line.split()
                
       
        # TODO need all snp locations in advance (2 passes)
        
