#!/usr/bin/bash

### capture sample name
sample=`basename $1`
sample=${sample%_L001_R1_001.fastq.gz}
sample=${sample%_L001_R2_001.fastq.gz}
echo 'sample: ' $sample
R1=${sample}_L001_R1_001.fastq.gz
R2=${sample}_L001_R2_001.fastq.gz

### capture analysis name
analysis=$2
echo 'analysis: ' $analysis

### BLAST pipeline
echo "ANALYSIS := $analysis" > 30-pipeline/blast/020-xml/analysis.make

make -C 30-pipeline/blast/010-fasta ${sample}_L001_R1_001.fasta
if [ -f "00-fastq/${sample}_L001_R2_001.fastq.gz" ]; then
    make -C 30-pipeline/blast/010-fasta ${sample}_L001_R2_001.fasta
else
    touch -r 30-pipeline/blast/010-fasta/${sample}_L001_R1_001.fasta 30-pipeline/blast/010-fasta/${sample}_L001_R2_001.fasta
fi
make -C 30-pipeline/blast/015-join $sample.fasta
make -C 30-pipeline/blast/020-xml $sample.xml
make -C 30-pipeline/blast/030-mut $sample.mut
