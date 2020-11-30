#!/bin/sh
# Runs the appropriate coverage summary perl script
#
# This must be run from the desired output summary directory
#
# $1 = segment type (WGS, Exon, CodingExon, or Intron)
# $2 = sample name, e.g.,SJRB001_D, SJMB002_D_102927016, SJRB001_D
# $3 = gender (M or F)
# $4 = presum directory
# $5 = chr2gene directory (ignored for WGS)

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
STYPE=$1  
SAMPLE=$2
GENDER=$3
PRESUM=$4
CHR2GENE=$5

case $STYPE in
  WG)           wg_summary.pl $GENDER $SAMPLE $PRESUM ;;
  Exon)         summary_by_segments.pl $GENDER $SAMPLE $PRESUM $CHR2GENE $STYPE ;;
  CodingExon)   summary_by_segments.pl $GENDER $SAMPLE $PRESUM $CHR2GENE $STYPE ;;
  Intron)       summary_by_segments.pl $GENDER $SAMPLE $PRESUM $CHR2GENE $STYPE ;;
  *)            echo "Invalid segment type $STYPE" >&2; exit 1 ;;
esac
