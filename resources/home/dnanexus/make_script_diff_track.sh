#!/usr/bin/env bash
#
# Creates coverage-diff commands for the given sample pair, and outputs
# those commands to stdout. 
#
# Params:
# $1 = Genome (Ex: GRCh37-lite)
# $2 = Case sample name (Ex: SJEWS001303_D1)
# $3 = Control sample name (Ex: SJEWS001303_G1)
# $4 = Input dir containing coverage wig files
# $5 = Output dir where results should be stored. This script will NOT create the output dir if it doesn't already exist.
# $6 = Coverage summary file (Ex: /nfs_exports/genomes/1/projects/WHOLEGENOME/ClinicalPilot/.BucketFinal/SJEWS/Coverage/Summary/SJEWS001303_D1.coverage.summary.table.txt)

# Show usage if insufficient number of params provided
if [ "$#" -eq 0 ]; then about.sh $0; exit 1; fi

# Set params
genome=$1
case_sample=$2
control_sample=$3
input_dir=$4
output_dir=$5
cov_summary_file=$6

# Get the mean coverage
case_smp_mean_cov=`cat $cov_summary_file | grep "wg.diagnosis" | cut -f8`
ctrl_smp_mean_cov=`cat $cov_summary_file | grep "wg.germline" | cut -f8`

# Make sure the mean coverage is not empty
if [ -z $case_smp_mean_cov ]
then
	echo "Error: No mean coverage for $case_sample. Exiting..."
	exit 1
fi
if [ -z $ctrl_smp_mean_cov ]
then
	echo "Error: No mean coverage for $control_sample. Exiting..."
	exit 1
fi

# Load the genome config
. import_config.sh genome $genome CHR_SIZES

# Create coverage-diff commands and print them to stdout
while read chr size
do
	echo "cd $output_dir; diff_track.pl ${input_dir}/${control_sample}_${chr}.wig ${input_dir}/${case_sample}_${chr}.wig $ctrl_smp_mean_cov $case_smp_mean_cov $chr"
done < $CHR_SIZES
