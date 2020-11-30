#!/bin/bash
# Runs Bambino's Ace2.SAMFinneyCoverage to compute coverage on a BAM and
# converts the results to bigwig (bw)
#
# Requires wigToBigWig to be in the PATH.
#
# $1 = input BAM file
# $2 = chromosome sizes file
# $3 = output bigwig file
# $4 = (optional) output gzipped wiggle file

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
INPUT=$1
CHR_SIZES=$2
OUTPUT=$3
OUTPUT_WIG=$4

# Write diagnostics
echo $0
hostname
date

# Make sure output directory exists
OUT_DIR=`dirname $OUTPUT`
if [ "$OUT_DIR" != "" -a ! -e "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi

# Create temp dir to hold intermediate files
SCRATCH=`mktemp -d` || ( echo "Could not get scratch dir" >&2 ; exit 1 )

# Intermediate gzipped wiggle
intermed=$SCRATCH/intermed.wig.gz

# Run Ace2.SAMFinneyCoverage
sfccmd="java.sh -Xmx8192m Ace2.SAMFinneyCoverage -bam $INPUT -wig -stdout"
gzcmd="gzip -1"
echo $sfccmd \| $gzcmd \> $intermed
$sfccmd | $gzcmd > $intermed
sfcec=${PIPESTATUS[0]}
if [ $sfcec -gt 0 ]
then
  echo "Ace2.SAMFinneyCoverage failed with exit code $sfcec" >&2
  exit 1
fi

# Debug
echo "Ace2.SAMFinneyCoverage complete.  Scratch dir contents:"
ls -l $SCRATCH

# Copy to output (optional)
if [ $OUTPUT_WIG ]
then
  echo "Copying gzipped wiggle file..."
  cp -v $intermed $OUTPUT_WIG
fi

# Convert to bigwig
echo "Converting to bigwig..."
L_OUTPUT=$SCRATCH/final.bw
cmd="wigToBigWig $intermed $CHR_SIZES $L_OUTPUT"
echo $cmd
if $cmd
then :
else
  echo "wigToBigWig failed" >&2
  exit 1
fi

echo "Conversion to bigwig complete.  Scratch dir contents:"
ls -l $SCRATCH

# Copy output
echo "Copying bigwig file..."
cp -v $L_OUTPUT $OUTPUT
