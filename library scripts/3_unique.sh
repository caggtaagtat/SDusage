#!/bin/bash
#PBS -l select=1:ncpus=8:mem=50gb
#PBS -l walltime=71:59:59
#PBS -r n

## Go to correct directory
GAPPATH=/path/to/library/bam/
cd $GAPPATH

## For every file in folder
for i in *.bam

## Do the following
do

## Get base name
i2=$(basename "${i%.bam}")

## Get only uniquely mapped reads for downstream analysis
/path/to/data/Apps/samtools-1.7_built/bin/samtools view -q 255 -b $i2.bam > $i2.q255.bam

## Do sorting and indexing of the new bam file
/path/to/data/Apps/samtools-1.7_built/bin/samtools sort -o $i2.q255.sorted.bam $i2.q255.bam 
/path/to/data/Apps/samtools-1.7_built/bin/samtools index $i2.q255.sorted.bam $i2.q255.sorted.bam.bai

## Finish the Loop Action
done

