#!/bin/bash

#############################################################################
# Script Name:  process_genome <accession> <id> <overwrite?> <log_file>      #
# Author:       Joel Hoefs .July 2023                                       #
# Description:  Downloads, prepares and gathers                             #
#               a single e.coli genome's variants into results/variants.    #
#               Designed only to be used through associated python script   #
#############################################################################

# special bash variable
SECONDS=0

cd ..

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
if [[ -f ${ACCESSION}_1.fastq ]] && [[ -f ${ACCESSION}_2.fastq ]] && [[ "$OVERWRITE" == "False" ]]; then
    echo "⚠️ fastq files ${ACCESSION}_x.fastq already exist, skipping .."
else
    echo "🔄 downloading pairwize 1 [${i}]"
    time wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${short}/${accession}/${accession}_1.fastq.gz
    echo "🔄 downloading pairwize 2 [${i}]"
    time wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${short}/${accession}/${accession}_2.fastq.gz
    echo "🔄 decompressing [${i}]"
    time yes n | gzip -d * 
fi
cd ..


### ALINGMENT ###
if [[ -f "results/${i}.bam" ]] && \ 
    [[ -f "results/${i}_sorted.bam" ]] && \ 
    [[ -f "results/${i}_sorted.bam.bai" ]] && \
    [[ "$OVERWRITE" == "False" ]]; then 
    echo "⚠️ alingmnet/sorted/indexed files for ${i}.bam already exists, skipping .."
else
    echo "🔄 aligning [${i}]"
    time bwa mem -M -t 2 \
        reference_data/ecoli_reference_k12 \
        raw_data/${accession}_2.fastq raw_data/${accession}_1.fastq \
        | samtools view -bS > results/${i}.bam;

    echo "🔄 sorting [${i}]"
    time samtools sort results/${i}.bam -O bam -o results/${i}_sorted.bam

    echo "🔄 indexing [${i}]"
    time samtools index results/${i}_sorted.bam
fi

rm raw_data/${accession}*
echo "⚠️ fastq raw data shoule be removed, looking at fastq: "
ls -l raw_data

### VARIANT CALLING ###
if [[ -f "results/${i}_calls.vcf.gz" ]] && [[ "$OVERWRITE" == "False" ]]; then 
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

echo "total time taken: $(($SECONDS / 60))min $(($SECONDS % 60))s"

rm raw_data/${accession}*
rm results/${i}*

printf "waiting 3 seconds before next genome .. \n";
sleep 3
