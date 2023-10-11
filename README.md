# AB-resistance-in-ecoli 🧬

Antibiotic Resistance (ABR) is a global burden to the future of modern medicine and has garnered multidisciplinary efforts to regulate, understand and surveil the phenomenon. Federal stewardship plans have encouraged bioinformaticians to research statistical alternatives to ABR testing and optimize workflows in genomic ABR studies. Following this endeavor, we evaluate machine learning (ML) methods against the most common unary track infection, Escherichia Coli (E. Coli) and 4 commonly proscribed antibiotics; Ciprofloxacin, Cefotaxime, Ceftazidime and Gentamicin. 




## Folder Structure


To see the main results and code I recommend referring to **model_iteration.ipynb** in "EDA and ML" or **scripts** in "Variant calling pipeline". 

- You can also view model_iteration.ipynb on google colab (github has a hard time loading it): https://colab.research.google.com/drive/1fEAUUCyM15sDIrPPUdDG0Eh6gCoD69Nt?usp=sharing




### EDA and ML
Contains expolatory data anaylises and model iterations + evaluations

- **model_iteration.ipynb - final SNP Ecoli preprocessing, model iteration & evaluation of ABR in 1000 ecoli genomes**
- ab_resistance_genes.ipynb - proof of concept of supervised machine learning in gene classification 
- Nigeria_Antibiotic_Susceptibility.ipynb - initial EDA task of ABR in Nigeria



### Variant Calling Pipeline
Contains the initial (unused) genome preprocessing pipeline for aws use.

- Deamon - AWS service running and diagnosis files
- Logs - logging for standard output 
- Reference_data - reference genome fasta and index files
- Raw_data - end point for raw downloaded fastq par end reads
- Results - bam, vcf and sorted, indexed files during preprocessing
    Variants - resultant snp and indel calls per genome
- **Scripts - python / bash scripts**
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

