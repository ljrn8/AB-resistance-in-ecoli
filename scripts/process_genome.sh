#!/bin/bash

#
# process_genome <accession> <id> <overwrite?> <log_file>     
# 
# Author:       Joel Hoefs .July 2023                                       
# Description:  Downloads, prepares and gathers                             
#               a single e.coli genome's variants into results/variants.    
#               Designed only to be used through associated python script   
#

# immediantly exit the script when any command fails
set -e

cd ..

SECONDS=0
accession=$1
i=$2
overwrite=$3
log_file=$4

short=${accession:0:6}
snp_file="results/variants/${i}_snps.txt"
indel_file="results/variants/${i}_indels.txt"

mkdir -p logs raw_data results results/variants;

if [[ -f "${snp_file}" ]] && [[ "$overwrite" = "False" ]]; then
    echo "⚠️ ${snp_file} already exists, skipping this genome .."
    exit 0
fi

printf "\n============================================\n"
printf "\tPROCESSING: ${i}, ${accession}\n"
printf "============================================\n"


### DOWNLOADING ###
cd raw_data
if [[ -f ${accession}_1.fastq ]] && [[ -f ${accession}_2.fastq ]] && [[ "$overwrite" == "False" ]]; then
    echo "⚠️ fastq files ${accession}_x.fastq already exist, skipping .."
else
    echo "🔄 downloading pairwize 1 [${i}]"
    time wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${short}/${accession}/${accession}_1.fastq.gz >/dev/null 2>&1
    echo "🔄 downloading pairwize 2 [${i}]"
    time wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${short}/${accession}/${accession}_2.fastq.gz >/dev/null 2>&1
    echo "🔄 decompressing [${i}]"
    time yes n | gzip -d ${accession}_1.fastq.gz ${accession}_2.fastq.gz
fi
cd ..


### ALINGMENT ###
if [[ -f "results/${i}.bam" ]] && [[ "$overwrite" == "False" ]]; then 
    echo "⚠️ alingmnet file for ${i}.bam already exists, skipping alingemnt .."
else
    echo "🔄 aligning [${i}]"
    time bwa mem -M -t 2 \
        reference_data/ecoli_reference_k12 \
        raw_data/${accession}_2.fastq raw_data/${accession}_1.fastq \
        | samtools view -bS > results/${i}.bam;
fi

# rm raw_data/${accession}*
echo "⚠️ fastq  files removed, raw_data folder: "
ls -l raw_data

echo "🔄 sorting [${i}]"
time samtools sort -m 100M results/${i}.bam -O bam -o results/${i}_sorted.bam

echo "🔄 indexing [${i}]"
time samtools index results/${i}_sorted.bam

### VARIANT CALLING ###
if [[ -f "results/${i}_calls.vcf.gz" ]] && [[ "$overwrite" == "False" ]]; then 
    echo "⚠️ variant calls already exist for ${i}_calls.vcf.gz, skipping .."
else
    echo "🔄 variant calling [${i}]"
    time bcftools mpileup --max-depth 500 -f reference_data/ecoli_reference_k12.fasta results/${i}_sorted.bam \
        | bcftools call -vm -Oz > results/${i}_calls.vcf.gz;
fi

echo "🔄 vcf cleaning [${i}]"
time bcftools view -Oz -e 'QUAL <= 20 || DP > 250 || MQBZ < -3 || RPBZ < -3 || RPBZ > 3 || SCBZ > 6' \
    results/${i}_calls.vcf.gz > results/${i}_filtered.vcf.gz
rm results/${i}_calls.vcf.gz

echo "🔄 collecting snps [${i}]"
time bcftools query "-i" 'TYPE="SNP"' -f '%POS %REF %ALT %QUAL\n' results/${i}_filtered.vcf.gz > "$snp_file"

echo "🔄 collecting indels [${i}]"
time bcftools query "-i" 'TYPE="INDEL"' -f '%POS %REF %ALT %QUAL\n' results/${i}_filtered.vcf.gz > "$indel_file"

wc -l $snp_file;
wc -l $indel_file;

rm results/${i}*

echo "total time taken: $(($SECONDS / 60))min $(($SECONDS % 60))s"

printf "waiting 3 seconds before next genome .. \n";
sleep 3
