import subprocess
import pandas as pd
import time 
from datetime import datetime
import os

# multi threading?
# 2118482 A C,G


def main(indexes=range(5)):
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
            bash process_genome.sh n {accession} {i} | tee -a {log_file};',
            shell=True
        )
        if exit_code != 0:
            print(f"non zero exit code from sub process {i}")
            return
        
        snp_file = f'../results/snps/{i}_snps.txt'
        if not os.path.exists(snp_file):
            print(f"snp file [{i}_snps] doesnt exists, assuming error")
            return
        
        print("waiting 4 seconds before next genome ... \n\n")
        time.sleep(4)


if __name__ =="__main__":
    print("this script downloads all genomes and may take few hours. \
        \nif you have cloned this repo and curious about its function, make sure to not ever run multiple instances of this program at once")
    x = input("Are you sure you want to proceed [Y/N]: ")
    if x in ['y', 'Y']:
        main()
    else:
        print("execution cancelled")

        