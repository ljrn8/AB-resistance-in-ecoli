#!/bin/bash

# TODO print more analytics
# TODO delete inside raw_data/

cd ..

ACCESSION=$1
i=$2
OVERWRITE=$3

SHORT=${ACCESSION:0:6}
SNP_FILE="results/snps/${i}_snps.txt"

mkdir -p logs raw_data results results/snps;

if [[ -f "${SNP_FILE}" ]] && [[ "$OVERWRITE" = "False" ]]; then
    echo "⚠️ ${SNP_FILE} already exists, skipping this genome .."
    exit 0
fi

printf "\n============================================\n"
printf "\tPROCESSING: ${i}, ${ACCESSION}\n"
printf "============================================\n"

### DOWNLOADING ###
cd raw_data;
if [[ -f ${ACCESSION}_1.fastq ]] && [[ -f ${ACCESSION}_2.fastq ]] && [[ "$OVERWRITE" = "False" ]]; then
    echo "⚠️ fastq files ${ACCESSION}_x.fastq already exist, skipping .."
else
    echo "🔄 downloading pairwize 1"
    wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SHORT}/${ACCESSION}/${ACCESSION}_1.fastq.gz
    echo "🔄 downloading pairwize 2"
    wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SHORT}/${ACCESSION}/${ACCESSION}_2.fastq.gz
    echo "🔄 decompressing"
    yes n | gzip -d * 
fi
cd ..;

### ALINGMENT ###
if [[ -f "results/${i}.bam" ]] && \ 
    [[ -f "results/${i}_sorted.bam" ]] && \ 
    [[ -f "results/${i}_sorted.bam.bai" ]] && \
    [[ "$OVERWRITE" = "False" ]]; then 
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
if [[ -f "results/${i}_calls.vcf.gz" ]] && [[ "$OVERWRITE" = "False" ]]; then 
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

printf "waiting 3 seconds before nect genome .. \n";
sleep 3