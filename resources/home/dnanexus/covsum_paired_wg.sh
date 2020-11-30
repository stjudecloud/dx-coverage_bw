#!/bin/sh
# Runs coverage summary paired processing (or the unpaired equivalent of the
# paired step)
#
# NOTE: THIS SCRIPT ASSUMES/ENFORCES CERTAIN CONVENTIONS.  You will specify a
# coverage directory and sample name, and the file locations will be:
# Input files: <dir>/presum/<sample>*
# Input BW file: <dir>/<subject of sample>/<sample>.bw
#
# Requires the pre-summary files *.Exon.txt and *.CodingExon.txt (and *.WG.txt 
# for WGS) should have been created at $COV_DIR/presum
#
# August 2012
#
# $1 case/tumor sample name, e.g. SJRB001_D, SJMB002_D_102927016
# $2 control/germline sample name, or "" or "-" for unpaired, e.g. SJRB001_G
# $3 gender, F for female and M for male, set to F if unknown
# $4 coverage file root directory e.g.
#    /nfs_exports/genomes/1/PCGP/HematopoieticMalignancies/SJTALL/NextGen/WholeGenome/Coverage

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
CASE_SAMPLE=$1
CONTROL_SAMPLE=$2
GENDER=$3
COV_DIR=$4
COV_CEILING=101

# Validate gender
if [ $GENDER != "F" -a $GENDER != "M" ]; then  echo "Invalid gender! Only 'F' and 'M' are allowed!" >&2; exit 1; fi

# Determine the output directory, create it if it doesn't exist, and change to
# it (the perl scripts write their output to the current working dir)
OUT_DIR=$COV_DIR/Summary
if [ $GENDER == "F" ]; then OUT_DIR=$OUT_DIR/F; fi
if [ ! -d $OUT_DIR ]; then mkdir -p $OUT_DIR; fi
cd $OUT_DIR

if [ $CONTROL_SAMPLE ]
then Rscript `dirname $0`/cov_hist.R $CONTROL_SAMPLE $CASE_SAMPLE $COV_CEILING
else Rscript `dirname $0`/unpaired_cov_hist.R  $CASE_SAMPLE $COV_CEILING
#then R CMD BATCH --no-save --no-restore "--args tumor=\"$CASE_SAMPLE\" normal=\"$CONTROL_SAMPLE\" covs=101" `dirname $0`/cov_hist.R
#else R CMD BATCH --no-save --no-restore "--args id=\"$CASE_SAMPLE\" covs=101" `dirname $0`/unpaired_cov_hist.R
fi
