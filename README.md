# AB-resistance-in-ecoli
Snp calling pipeline for predicting Antibiotic resistance in E.coli from SNP's. Below is a description of the folder strucuture. 
To see the main results and code I recomend referring to "EDA and ML"/TODO or "Variant calling pipeline"/Scripts

## TODO
- add final notebook and change to links + add diagrams to README


## EDA and ML
Contains expolatory data anaylises and model iterations + evaluations

- **TODO - final SNP Ecoli preprocessing, model iteration & evaluation of ABR in 1000 ecoli genomes**
- ab_resistance_genes.ipynb - proof of concept of supervised machine learning in gene classification 
- Nigeria_Antibiotic_Susceptibility.ipynb - initial EDA task of ABR in nigeria


## variant calling pipeline
Contains the initial (unused) genome preprocessing pipeline for aws use.

- Deamon - AWS service running and diagnosis files
- Logs - logging for standard output ()
- Reference_data - reference genome fasta and index files
- Raw_data - end point for raw downlaoded fastq par end reads
- Results - bam, vcf and sorted, indexed files during preprocessing
    Variants - resultant snp and indel calls per genome
- **Scripts - python / bash scripts**
    (preprocess_variants.py is not used)
- Figures (untracked)
    Output diagrams 
- EDA (untracked) preliminary data analysis




---
### Dependancies
- python3
    - pandas
- bwa
- bash
- samtools
- bcftools
- wget

