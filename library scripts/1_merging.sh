#!/bin/bash
#PBS -l select=1:ncpus=4:mem=20gb
#PBS -l walltime=70:59:59
#PBS -r n

## Change to correct working directory
cd /path/to/library/

## Merge forward and reverse reads for data of both plasmid and SD usage samples
/path/to/data/Apps/bbmap/bbmerge.sh in1=Downstream_SD_usage_R1.fastq in2=Downstream_SD_usage_R2.fastq out=Downstream_SD_usage_merged.fq outu=Downstream_SD_usage_unmerged.fq

/path/to/data/Apps/bbmap/bbmerge.sh in1=Plasmid_R1.fastq in2=Plasmid_R2.fastq out=Plasmid_merged.fq outu=Plasmid_unmerged.fq
