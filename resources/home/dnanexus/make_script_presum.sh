#!/bin/bash
# Runs post-coverage steps to prepare for coverage summary.
#
# NOTE: THIS SCRIPT ASSUMES/ENFORCES CERTAIN CONVENTIONS.  You will specify a
# coverage directory and subject/sample names, and the file locations will be:
# Input BW file: <dir>/<subdir>/<sample>.bw
# Output files: <dir>/presum/<sample>*
#
# Gang Wu, August 7th, 2012
#
# Example sh presum.sh hg19 /nfs_exports/genomes/1/PCGP/.allTumor/SJMB/NextGen/WholeGenome/Coverage SJMB023_D SJMB023
#
# Modified to use config.sh to parse the CHR2GENE annotation path by genome
# version for future extension to different species
#
# Required: bigwig coverage files under $COV_DIR/$DONOR/
# Note: some long chromosomes (such as hg19 chr1, chr2 and chr3) may require up
# to 20 Gb memory to process, others should be done <15Gb
#
# $1 = genome version, case sensitive (e.g. hg19, hg18, mm9)
# $2 = coverage directory (SEE NOTE ON CONVENTIONS ABOVE), e.g.
#      /nfs_exports/genomes/1/PCGP/.allTumor/SJMB/NextGen/WholeGenome/Coverage
# $3 = subject ID, e.g., SJMB040
# $4 = sample name (pid only, not full barcode), e.g., SJMB040_D or SJMB040_G
# $5,... = Segment types, e.g., WG CodingExon Exon for WGS

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
GENOME=$1
COV_DIR=$2
SUBJECT=$3
SAMPLE=$4
shift 4
SEGMENT_TYPES="$@"

# Find the input bw file
BW=$COV_DIR/$SUBJECT/$SAMPLE.bw

# Load genome config
. import_config.sh genome $GENOME

# Create the output directory if it does not exist
if [ ! -d $COV_DIR/presum ]; then mkdir -p $COV_DIR/presum; fi

# Loop through all canonical chromosomes (those in CHR_SIZES), stripping off
# the "chr" prefix
for chr in `cat $CHR_SIZES | cut -f 1 | sed 's/^chr//'`;do
  #dont do chrM and chrMT
  if [ "$chr" == "M" -o "$chr" == "MT" ]; then continue; fi
  
  # Loop over the segment types
  for STYPE in $SEGMENT_TYPES
  do
    # Build the command
    cmd="bw_presum.pl -chr $chr -cov_file $BW -annot_dir $CHR2GENE -chr_size $CHR_SIZES -out_dir $COV_DIR/presum -sample $SAMPLE -seg_type $STYPE"
    # Write it out
    echo $cmd
  done
done
