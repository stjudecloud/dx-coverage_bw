#!/bin/bash
# Writes commands to run bw_presum.pl
#
# NOTE: THIS SCRIPT REPLACES make_script_presum.pl, WHICH ENFORCED CERTAIN
# DIRECTORY STRUCTURE CONVENTIONS.
#
# Gang Wu, August 7th, 2012
#
# Example sh presum.sh hg19 /nfs_exports/genomes/1/PCGP/.allTumor/SJMB/NextGen/WholeGenome/Coverage SJMB023_D SJMB023
#
# Modified to use config.sh to parse the CHR2GENE annotation path by genome
# version for future extension to different species
#
# Required: bigwig coverage files
# Note: some long chromosomes (such as hg19 chr1, chr2 and chr3) may require up
# to 20 Gb memory to process, others should be done <15Gb
#
# $1 = genome version, case sensitive (e.g. GRCh37-lite, MGSCv37)
# $2 = input directory
# $3 = output directory
# $4 = sample name (pid only, not full barcode), e.g., SJMB040_D or SJMB040_G
# $5,... = Segment types, e.g., WG CodingExon Exon for WGS

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
GENOME=$1
IN_DIR=$2
OUT_DIR=$3
SAMPLE=$4
shift 4
SEGMENT_TYPES="$@"

# Find the input bw file
BW=$IN_DIR/$SAMPLE.bw

# Load genome config
. import_config.sh genome $GENOME

# Create the output directory if it does not exist
if [ ! -d $OUT_DIR ]; then mkdir -p $OUT_DIR; fi

# Loop through all canonical chromosomes (those in CHR_SIZES), stripping off
# the "chr" prefix
for chr in `cat $CHR_SIZES | cut -f 1 | sed 's/^chr//'`;do
  #dont do chrM and chrMT
  if [ "$chr" == "M" -o "$chr" == "MT" ]; then continue; fi
  
  # Loop over the segment types
  for STYPE in $SEGMENT_TYPES
  do
    # Build the command
    cmd="bw_presum.pl -chr $chr -cov_file $BW -annot_dir $CHR2GENE -chr_size $CHR_SIZES -out_dir $OUT_DIR -sample $SAMPLE -seg_type $STYPE"
    # Write it out
    echo $cmd
  done
done
