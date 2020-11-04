#!/usr/bin/bash

echo Start on $(date --iso-8601='seconds')

[ ! -d rawdata ] && mkdir rawdata
cd rawdata

echo
echo 1 Downloading and recompressing data

if [ ! -f GCF_000005845.2_ASM584v2_genomic.fna ]
then
    [ ! -f GCF_000005845.2_ASM584v2_genomic.fna.gz ] && wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
    
    gzip -d GCF_000005845.2_ASM584v2_genomic.fna.gz
fi
if [ ! -f GCF_000005845.2_ASM584v2_genomic.gff ]
then    
    [ ! -f GCF_000005845.2_ASM584v2_genomic.gff.gz ] && wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz

    gzip -d GCF_000005845.2_ASM584v2_genomic.gff.gz
fi

if [ ! -f amp_res_1.fastq.gz ]
then
    if [ ! -f amp_res_1.fastq ]
    then 
        [ ! -f amp_res_1.fastq.zip ] && wget http://public.dobzhanskycenter.ru/mrayko/amp_res_1.fastq.zip
        unzip amp_res_1.fastq.zip && rm amp_res_1.fastq.zip 
    fi
    gzip --verbose amp_res_1.fastq
fi
if [ ! -f amp_res_2.fastq.gz ]
then
    if [ ! -f amp_res_2.fastq ]
    then 
        [ ! -f amp_res_2.fastq.zip ] && wget http://public.dobzhanskycenter.ru/mrayko/amp_res_2.fastq.zip
        unzip amp_res_2.fastq.zip && rm amp_res_2.fastq.zip 
    fi
    gzip --verbose amp_res_2.fastq
fi

echo
echo 2 Inspect raw data
echo $(($(zcat amp_res_1.fastq.gz | wc -l) / 4)) reads

echo
echo 3 Basic statistic

# http://drive5.com/usearch/manual/quality_score.html

fastqc -o . amp_res_?.fastq.gz

# http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/

echo
echo 4 Filtering the reads

# http://www.usadellab.org/cms/?page=trimmomatic

TrimmomaticPE -trimlog amp_res_trimming.log -phred33 amp_res_1.fastq.gz amp_res_2.fastq.gz amp_res_paired_output_1.fq.gz amp_res_unpaired_output_1.fq.gz amp_res_paired_output_2.fq.gz amp_res_unpaired_output_2.fq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20

# Input Read Pairs: 455876 Both Surviving: 446259 (97.89%) Forward Only Surviving: 9216 (2.02%) Reverse Only Surviving: 273 (0.06%) Dropped: 128 (0.03%)

echo $(($(zcat amp_res_paired_output_1.fq.gz | wc -l) / 4)) reads

echo
echo Basic statistics of filtered data

fastqc -o . amp_res_paired_output_?.fq.gz

echo
echo 5 Aligning sequences to reference

echo 5.1 Index the reference
bwa index GCF_000005845.2_ASM584v2_genomic.fna
echo 5.2 Align reads
bwa mem GCF_000005845.2_ASM584v2_genomic.fna amp_res_paired_output_1.fq.gz amp_res_paired_output_2.fq.gz > alignment.sam

echo 5.3 Compress and sort sam to bam
samtools view -S -b alignment.sam > alignment.bam
echo Some basic statistics
samtools flagstat alignment.bam
echo 5.4 Sort bam by sequence coordinate on reference
samtools sort alignment.bam -o alignment_sorted.bam
samtools index alignment_sorted.bam

echo
echo 6 Variant calling
echo Pile up the reads
samtools mpileup -f GCF_000005845.2_ASM584v2_genomic.fna alignment_sorted.bam > my.mpileup

echo Call actual variants
java -jar ~/software/VarScan.v2.3.9.jar mpileup2snp my.mpileup --min-var-freq 0.5 --variants --output-vcf 1 > VarScan_results.vcf

# Reference genome should be first registered in snpEff.config and build with
# java -jar snpEff.jar build -gff3 ecoli
# echo
# echo Annotate variants
# java -jar snpEff.jar eff ecoli VarScan_results.vcf > ann.txt

echo
echo End on $(date --iso-8601='seconds')
