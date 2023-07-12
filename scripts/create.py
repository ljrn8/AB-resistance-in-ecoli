import subprocess
import pandas as pd
import time 
from datetime import datetime
import os, sys

# TODO can be done with maximum 3 threads
# 2118482 A C,G

def main(
    overwrite=False,
    indexes=range(5)
):
    now = datetime.now()
    now_str = now.strftime("%d-%m-%Y_%H:%M:%S")
    df = pd.read_csv("../accessions.csv")
    log_file = f"../logs/{now_str}_log.txt"
    
    for i in indexes:
        accession = df["Lane.accession"].iloc[i]

        # genome setup script
        exit_code = subprocess.call(f' \
            touch {log_file}; \
            bash process_genome.sh {accession} {i} {overwrite} | tee -a {log_file};',
            shell=True  
        )
        if exit_code != 0:
            print(f"non zero exit code from sub process {i}")
            return 
        
        snp_file = f'../results/variants/{i}_snps.txt'
        if not os.path.exists(snp_file):
            print(f"snp file [{i}_snps] doesnt exists, assuming error")
            return
        

from argparse import ArgumentParser

if __name__ =="__main__":

    # index flags
    parser = ArgumentParser()
    parser.add_argument("-s", "--start", default=0, type=int, help="starting index")
    parser.add_argument("-e", "--ending", default=5, type=int, help="ending index (exclusive)")
    parser.add_argument("-o","--overwrite",action="store_true", help="overwrite all files") 
    parser.add_argument("-y","--yes",action="store_true", help="skip yessing") 

    args = vars(parser.parse_args()) 
    start, end = args["start"], args["ending"]
    if start > end:
        raise ValueError("irrational index range")

    overwrite = args["overwrite"]
    if overwrite:
        print("⚠️ overwriting files '-o'")

    index_range = range(start, end)
    print(f'(using index range {index_range}')

    print("\nthis script downloads all genomes and may take few hours. \
        \nif you have cloned this repo and curious about its function, make sure to not ever run multiple instances of this program at once")
    
    if not args["yes"]: 
        x = input("Are you sure you want to proceed [Y/N]: ")
        if x in ['y', 'Y']:
            main(overwrite=overwrite, indexes=index_range)
    else:
        main(overwrite=overwrite, indexes=index_range)
    
    print("execution finished")

        
