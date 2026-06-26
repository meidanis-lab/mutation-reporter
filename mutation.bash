#!/usr/bin/bash

### capture sample name
fastq=`basename $1`
if [[ "$fastq" =~ R1_001.fastq ]]; then
    R1fastq=$fastq
    R1fasta=${R1fastq/fastq.gz/fasta}
    R2fastq=${R1fastq/R1_001.fastq/R2_001.fastq}
    R2fasta=${R2fastq/fastq.gz/fasta}
else
    R2fastq=$fastq
    R2fasta=${R2fastq/fastq.gz/fasta}
    R1fastq=${R2fastq/R2_001.fastq/R1_001.fastq}
    R1fasta=${R1fastq/fastq.gz/fasta}
fi
sample=${R1fastq%.fastq.gz}
sample=${sample%_R1_001}
sample=${sample%_R2_001}
sample=${sample%_L001}
echo 'sample: ' $sample

### capture analysis name
analysis=$2
echo 'analysis: ' $analysis

### BLAST pipeline
echo "ANALYSIS := $analysis" > 30-pipeline/blast/020-xml/analysis.make

cat 30-pipeline/blast/010-fasta/Makefile
make -C 30-pipeline/blast/010-fasta $R1fasta
echo 'made R1 fasta'
if [ -f "00-fastq/$R2fastq" ]; then
    make -C 30-pipeline/blast/010-fasta $R2fasta
else
    touch -r 30-pipeline/blast/010-fasta/$R1fasta 30-pipeline/blast/010-fasta/$R2fasta
fi
echo 'made R2 fasta'
make -C 30-pipeline/blast/015-join $sample.fasta
echo 'made join'
make -C 30-pipeline/blast/020-xml $sample.xml
echo 'made blast'
make -C 30-pipeline/blast/030-mut $sample.mut
echo 'made mut'
