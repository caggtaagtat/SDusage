#!/bin/bash
#PBS -l select=1:ncpus=8:mem=40gb
#PBS -l walltime=70:59:59
#PBS -r n

## Go to directory with annotations 
cd /path/to/data/annotations/human/

## Fuse the plasmid.fasta file with the fasta file of the human reference genome GRCh38 
#awk 1 plasmid.fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa > hg38_plasmid_fusion.fasta

## Do the same for the gtf file
#cat /path/to/data/annotations/human/plasmid.gtf Homo_sapiens.GRCh38.100.gtf > hg38_plasmid_fusion_annotation.gtf

## Execute STAR to generate a STAR index of the combined references
/path/to/data/Apps/STAR/STAR-2.5.4b/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /path/to/data/annotations/human/star_index_hg38_plasmid_r150/ --genomeFastaFiles hg38_plasmid_fusion.fasta --sjdbGTFfile hg38_plasmid_fusion_annotation.gtf --sjdbOverhang 149


## Change directory to fastq file directory
GAPPATH=/path/to/library/
cd $GAPPATH

## Unzip zipped files
gunzip *.gz

## For every samples whose file name is listed seperated by a new line in the text file "xaa"
while read SAMPLE; do

## Get the sole file name, without the path and file ending
FILEBASE=$(basename ${SAMPLE%"_merged.fq"})
  
## Make new directory for the sample
mkdir $GAPPATH/$FILEBASE.STAR

## Enter the new directory
cd $GAPPATH/$FILEBASE.STAR

## Align with STAR aligner
/path/to/data/Apps/STAR/STAR-2.5.4b/bin/Linux_x86_64/STAR --outFilterType BySJout --outFilterMismatchNmax 5 --outFilterMismatchNoverLmax 0.04 --alignEndsType EndToEnd --runThreadN 8 --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 4 --alignIntronMax 300000 --alignSJoverhangMin 8 --genomeDir /path/to/data/annotations/human/star_index_hg38_plasmid_r150/ --sjdbOverhang 149 --sjdbGTFfile /path/to/data/annotations/human/hg38_plasmid_fusion_annotation.gtf --outFileNamePrefix $GAPPATH/$FILEBASE.STAR/ --readFilesIn $SAMPLE > STARaligning.log 

done < $GAPPATH/xaa
