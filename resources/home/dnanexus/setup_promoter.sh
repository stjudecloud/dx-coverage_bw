#!/bin/sh
#
# Sets up scripts to compute coverage in promoter regions.
#
# $1 = Coverage root directory (Ex: /nfs_exports/genomes/1/projects/WGS/ClinicalPilot/.BucketFinal)
# $2 = File containing paths to BAM files
# $3 = Genome (Ex: hg19, hg18, mm9)
# $4 = 1 (for WGS) or 0 (for others)
#
# Author: Gang Wu
#
# TODO Use config system for genome-specific values
# TODO Use commands file/bsub_array_for_cmdfile.sh system

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Params
COVROOTDIR=$1
LIST=$2
GENOME=$3
WGS=$4
SCRIPT=presum_prep.sh

# Compute the coverage for each BAM file
for BAM in `cat $LIST`; do
  BAM_NAME=`basename $BAM`
  DISEASE=`sn_parse.sh $BAM_NAME | cut -d " " -f2`
  COVDIR=$COVROOTDIR/SJ$DISEASE/Coverage
  mkdir -p $COVDIR/presum;
  SAMPLE=`basename $BAM | cut -f1 -d "-"`
  DONOR=`sn_parse.sh $SAMPLE | sed 's/ /|/' | cut -f1 -d"|"`;
  $SCRIPT $GENOME $COVDIR ${SAMPLE} ${DONOR}

  ls ${SAMPLE}*.Promoter | perl -ne '`bsub -P PCGP -M 2000 -R "rusage[mem=2000]" <$_`'
done
