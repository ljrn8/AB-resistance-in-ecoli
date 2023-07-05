
# TODO bcf is for binary call format, consider vcf --gzvfc in gh
# TODO print analytics
# TODO delete inside raw_data/

ACCESSION=$1
SHORT=${ACCESSION:0:6}
i=$2

echo "======================================="
echo "PROCESSING: ${i}, ${ACCESSION}"

echo "downloading pairwize 1"
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SHORT}/${ACCESSION}/${ACCESSION}_1.fastq.gz -P ../raw_data/
echo "downloading pairwize 2"
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SHORT}/${ACCESSION}/${ACCESSION}_2.fastq.gz -P ../raw_data/
gzip -d raw_data/*


echo "STARTING ALINGMENT ..."
bwa mem -M -t 2 -R \
    reference_data/ecoli_reference_k12 \
    raw_data/${ACCESSION}_2.fastq raw_data/${ACCESSION}_1.fastq \
    | samtools view -Sb - > results/${i}.bam 


echo "SORTING ..."
samtools sort ../results/${i}.bam -O bam -o ../results/${i}_sorted.bam

echo "INDEXING ..."
samtools index ../results/${i}_sorted.bam

echo "VARIANT CALLING ..."
bcftools mpileup --max-depth 500 -f \
    reference_data/ecoli_reference_k12.fasta results/${i}_sorted.bam 
    | bcftools call -vm -Oz > results/${i}_calls.vcf.gz

echo "VCF CLEANING ..."
bcftools view -e 'QUAL <= 20 || DP > 250 || MQBZ < -3 || RPBZ < -3 || RPBZ > 3 || || SCBZ > 6' ${i}_calls.vcf > filtered.vcf

echo "COLLECTING SNPS ..."
bcftools query -i 'TYPE="SNP"' -f '%POS %REF %ALT\n' ${i}_calls.vcf > ${i}_snps

echo "======================================="