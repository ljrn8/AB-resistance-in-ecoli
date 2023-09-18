# AB-resistance-in-ecoli
Snp calling pipeline for predicting Antibiotic resistance in E.coli from SNP's 

### Folder structure
- Deamon - AWS service running and diagnosis files
- Logs - logging for standard output ()
- Reference_data - reference genome fasta and index files
- Raw_data - end point for raw downlaoded fastq par end reads
- Results - bam, vcf and sorted, indexed files during preprocessing
    Variants - resultant snp and indel calls per genome
- Scripts - python / bash scripts
    (preprocess_variants.py is not used)
- Figures (untracked)
    Output diagrams 
- EDA (untracked)


### Dependancies
- python3
    - pandas
- bwa
- bash
- samtools
- bcftools
- wget

