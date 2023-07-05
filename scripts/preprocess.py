
from Bio import SeqIO, Entrez

def download_genome(accession):
    
    
    handle = Entrez.efetch(db="nucleotide", id=accession)
    records = handle.readlines()
    
    
    
    for record in records:
        fields = record.strip().split(",")
        run_accession = fields[0]
        fastq_ftp = fields[9]
        ftp_links = fastq_ftp.split(";")


        for link in ftp_links:
            print(link)
            """ if link.endswith("_1.fastq.gz"):
                SeqIO.write(SeqIO.parse(link, "fastq"), run_accession + "_1.fastq", "fastq")
            elif link.endswith("_2.fastq.gz"):
                SeqIO.write(SeqIO.parse(link, "fastq"), run_accession + "_2.fastq", "fastq")
            """


def main() -> None:
    download_genome("ERS356929")

if __name__ =="__main__":
    # add warning about long run duration and overriding data
    
    print("this script downloads all 2000 genomes and should be ran sparingly")
    x = input("Are you sure you want to proceed [Y/N]: ")
    if x in ['y', 'Y']:
        main()
    
    print("execution cancelled")

        
