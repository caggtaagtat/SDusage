#!/bin/bash
#PBS -l select=1:ncpus=2:mem=10gb
#PBS -l walltime=71:59:59
#PBS -r n

## Go to correct directory
GAPPATH=/path/to/library/bam/
cd $GAPPATH

## For every file in the given path, showing a ".q255.sorted.bam" at the end of the filename
find $GAPPATH -name "*.q255.sorted.bam" | while read i

## Do the following
do

## Get the sole file name, without the path and ".q255.sorted.bam" file ending
FILEBASE=$(basename "${i%.q255.sorted.bam}")

## Count how many reads aligned to what "chromosome"
/path/to/data/Apps/samtools-1.7_built/bin/samtools idxstats $i > $FILEBASE.stats

## Finish the Loop Action
done


