#!/bin/bash

# TODO print more analytics
# TODO delete inside raw_data/

cd ..


# TODO override flag
ACCESSION=$2
i=$3

SHORT=${ACCESSION:0:6}
OVERRIDE=0
SNP_FILE="results/snps/${i}_snps.txt"

mkdir -p logs raw_data results results/snps;

while getopts "f" flag; do
    case "$flag" in
        f) 
            echo "⚠️ overriding ${SNP_FILE}"
            OVERRIDE=1
            ;;
    esac
done

if [ -f $SNP_FILE ] && [ $OVERRIDE -eq 0 ]; then
    echo "${SNP_FILE} already exists, run with -f to override snp files"
    exit 0
fi

printf "\n============================================\n"
printf "\tPROCESSING: ${i}, ${ACCESSION}\n"
printf "============================================\n"

### DOWNLOADING ###
echo "🔄 downloading pairwize 1"
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SHORT}/${ACCESSION}/${ACCESSION}_1.fastq.gz -P raw_data/
echo "🔄 downloading pairwize 2"
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SHORT}/${ACCESSION}/${ACCESSION}_2.fastq.gz -P raw_data/
echo "🔄 decompressing"
yes n | gzip -d raw_data/* 


### ALINGMENT ###
if [ -f "results/${i}.bam" ] && [ -f "results/${i}_sorted.bam" ] &&  [ -f "results/${i}_sorted.bam.bai" ]; then 
    echo "⚠️ alingmnet/sorted/indexed files for ${i}.bam already exists, skipping .."
else
    echo "🔄 aligning"
    bwa mem -M -t 2 \
        reference_data/ecoli_reference_k12 \
        raw_data/${ACCESSION}_2.fastq raw_data/${ACCESSION}_1.fastq \
        | samtools view -bS > results/${i}.bam;

    echo "🔄 sorting"
    samtools sort results/${i}.bam -O bam -o results/${i}_sorted.bam

    echo "🔄 indexing"
    samtools index results/${i}_sorted.bam
fi

### VARIANT CALLING ###
if [ -f "results/${i}_calls.vcf.gz" ]; then 
    echo "⚠️ variant calls already exist for ${i}_calls.vcf.gz, skipping .."
else
    echo "🔄 variant calling"
    bcftools mpileup --max-depth 500 -f reference_data/ecoli_reference_k12.fasta results/${i}_sorted.bam \
        | bcftools call -vm -Oz > results/${i}_calls.vcf.gz;
fi

echo "🔄 vcf cleaning"
bcftools view -Oz -e 'QUAL <= 20 || DP > 250 || MQBZ < -3 || RPBZ < -3 || RPBZ > 3 || SCBZ > 6' \
    results/${i}_calls.vcf.gz > results/${i}_filtered.vcf.gz

echo "🔄 collecting snps"
bcftools query "-i" 'TYPE="SNP"' -f '%POS %REF %ALT\n' results/${i}_filtered.vcf.gz > "$SNP_FILE"

echo "number of snps: "
wc -l $SNP_FILE;

printf "\n";