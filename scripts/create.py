import subprocess
import pandas as pd
import time 
from datetime import datetime
import os, sys

# TODO can be done with 3 threads
# look at AWS genomics & azure
# multi threading
# 2118482 A C,G


def main(
    overwrite=False,
    indexes=range(5)
):
    now = datetime.now()
    now_str = now.strftime("%d-%m-%Y_%H:%M:%S")
    df = pd.read_csv("../accessions.csv")
    log_file = f"../logs/{now_str}_log.txt"
    
    
    # just do 5 for now
    for i in indexes:
        accession = df["Lane.accession"].iloc[i]
        # TODO phenotypes
        # TODO time

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
        


if __name__ =="__main__":

    overwrite = "-o" in sys.argv
    if overwrite:
        print("⚠️ overwriting files '-o'")

    print("this script downloads all genomes and may take few hours. \
        \nif you have cloned this repo and curious about its function, be aware this is not for general use.")
    x = input("Are you sure you want to proceed [Y/N]: ")
    if x in ['y', 'Y']:
        main(overwrite=overwrite)
    else:
        print("execution cancelled")

        