#!/bin/sh
# Runs coverage summary
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
# $1 sample name, e.g.,SJRB001_D, SJMB002_D_102927016, SJRB001_D
# $2 gender, F for female and M for male, set to F if unknown
# $3 coverage file root directory e.g.
#    /nfs_exports/genomes/1/PCGP/HematopoieticMalignancies/SJTALL/NextGen/WholeGenome/Coverage
# $4 genome version, case sensitive (e.g. hg19, hg18 or mm9)
# $5,... = Segment types, e.g., "WG CodingExon Exon" for WGS

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
SAMPLE=$1  
GENDER=$2
COV_DIR=$3
# Doesn't depend on environment genome
CHR_SIZES=$4
shift 4
SEGMENT_TYPES="$@"

# Validate gender
if [ $GENDER != "F" -a $GENDER != "M" ]; then  echo "Invalid gender! Only 'F' and 'M' are allowed!" >&2; exit 1; fi

# Load genome config
#. import_config.sh genome $GENOME

CHR_NUM=`grep -v M $CHR_SIZES|wc -l`;
echo "# of chromosomes except chrM or chrMT: $CHR_NUM"

# Determine the presum directory
PRESUM=$COV_DIR/presum

# Determine the output directory, create it if it doesn't exist, and change to
# it (the perl scripts write their output to the current working dir)
OUT_DIR=$COV_DIR/Summary
if [ $GENDER == "F" ]; then OUT_DIR=$OUT_DIR/F; fi
if [ ! -d $OUT_DIR ]; then mkdir -p $OUT_DIR; fi
cd $OUT_DIR

# Loop over the segment types
for STYPE in $SEGMENT_TYPES
do
  # Count non-empty files
  FILENUM=`grep -E -m 1 'Coverage|Segment' $PRESUM/$SAMPLE*.$STYPE.txt|wc -l`
  if [ $FILENUM -lt 1 ]; then echo "There are no ${SAMPLE}*.$STYPE.txt under $PRESUM/" >&2; exit 1; fi
  MOD=`expr $FILENUM % $CHR_NUM`
  if [ $MOD -gt 0 ]; then echo "The number of ${SAMPLE}*.$STYPE.txt is not ${CHR_NUM}x!\n"; exit 1; fi
  
  # Run the summary script
  if ! run_covsum_pl.sh $STYPE $SAMPLE $GENDER $PRESUM $CHR2GENE
  then
    echo "run_covsum_pl.sh failed"
    exit 1
  fi

  # Run the cumulative coverage R script
  lcstype=`echo $STYPE | tr '[:upper:]' '[:lower:]'`
  if ! Rscript `dirname $0`/get_cumulative_coverage.R $SAMPLE.$lcstype
  then 
    echo "get_cumulative_coverage.R failed"
    exit 1
  fi

  # Run the boxplots R script
  cov_means_file=$SAMPLE.$lcstype.coverage.means.with.GC.txt
  if [ -e "$cov_means_file" ]
  then
    Rscript `dirname $0`/cov_boxplots.R $SAMPLE $lcstype
    exit_code=$?
    if [ "$exit_code" -ne 0 ]
    then 
      echo "cov_boxplots.R failed"
      exit 1
    fi
  fi
done

echo "covsum done"
date
