#!/usr/bin/env bash
#
# Creates commands to convert bigwig files to wig files for each chromosome
# and outputs those commands to stdout.
#
# Params:
# $1 = Genome (Ex: GRCh37-lite)
# $2 = Sample name (Ex: SJMB009_D)
# $3 = Input dir containing bigwig files
# $4 = Output directory where wig files will be stored. Output dir will NOT be created if it doesn't already exist

# Show usage if insufficient number of params provided
if [ "$#" -eq 0 ]; then about.sh $0; exit 1; fi

# Set Params
genome=$1
sample=$2
input_dir=$3
output_dir=$4

# Load the genome config
. import_config.sh genome $genome CHR_SIZES

# Create commands for each chr
while read chr size
do
	echo "bigWigToWig -chrom=$chr -start=0 -end=${size} $input_dir/${sample}.bw $output_dir/${sample}_${chr}.wig"
done < $CHR_SIZES
