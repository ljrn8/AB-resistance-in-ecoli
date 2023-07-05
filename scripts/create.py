

import subprocess
import pandas as pd
import time 
from datetime import datetime


def main():
    now = datetime.now()
    now_str = now.strftime("%d-%m-%Y_%H:%M:%S")
    df = pd.read_csv("../accessions.csv")
    
    # just do 5 for now
    for i in range(5):
        accession = df["Lane.accession"].iloc[i]
        # TODO phenotypes
        # TODO time
        
        subprocess.call(f' \
            touch ../logs/{now_str}_log.txt; \
            time bash process_genome.sh {accession} {i} | tee -a ../logs/{now_str}_log.txt;',
            shell=True
        )
        time.sleep(4)


if __name__ =="__main__":
    print("this script downloads all genomes and may override data \
        if you have cloned this repo and curious about its function, be aware of its many dependancies")
    x = input("Are you sure you want to proceed [Y/N]: ")
    if x in ['y', 'Y']:
        main()
    else:
        print("execution cancelled")

        