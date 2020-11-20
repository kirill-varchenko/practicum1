#!/usr/bin/bash

echo Start on $(date --iso-8601='seconds')

cd rawdata

bwa index sequence.fasta

files=( SRR1705851 SRR1705858 SRR1705859 SRR1705860 )
for i in "${files[@]}"
do
	echo $i
    bwa mem sequence.fasta $i.fastq.gz | samtools view -S -b - | samtools sort - -o ${i}_sorted.bam
    samtools index ${i}_sorted.bam
    # samtools view -f 4 SRR1705851_sorted.bam | samtools fasta - | grep -c '>'
    samtools mpileup -d 0 -f sequence.fasta ${i}_sorted.bam > $i.mpileup
    java -jar ~/software/VarScan.v2.3.9.jar mpileup2snp $i.mpileup --min-var-freq 0.95 --variants --output-vcf 1 > ${i}_95.vcf
    java -jar ~/software/VarScan.v2.3.9.jar mpileup2snp $i.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > ${i}_01.vcf
    cat ${i}_95.vcf | awk 'NR>24 {split($10, a, ":"); print $1, $2, $4, $5, a[7]}' > ${i}_95.csv
    cat ${i}_01.vcf | awk 'NR>24 {split($10, a, ":"); print $1, $2, $4, $5, a[7]}' > ${i}_01.csv
done


echo
echo End on $(date --iso-8601='seconds')
