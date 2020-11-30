#!/bin/bash
#
# Computes coverage in WG, exons, coding exons, introns, and promoters.
# 
# $1 = Genome (hg19, hg18, mm9)
# $2 = Path to coverage directory
# $3 = Sample name (Ex: SJETV092_D)
# $4 = Donor (Ex: SJETV092)
# 
# Author: Gang Wu, June 1st, 2011
#
# TODO Use config system for genome-specific values
# TODO Use commands file/bsub_array_for_cmdfile.sh system

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Params
GENOME=$1
COV_DIR=$2
SAMPLE=$3
DONOR=$4
BW=$COV_DIR/$DONOR/$SAMPLE.bw

# Get chromosome size file and annotation file 
if [ $GENOME == "hg19" ]; then
  CHRSIZE=/nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/GRCh37.sizes
  ANNOT=/nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/chr2gene
elif [ $GENOME == "hg18" ]; then
  CHRSIZE=/nfs_exports/apps/gnu-apps/NextGen/nextgensupport/hg18/hg18.sizes
  ANNOT=/nfs_exports/apps/gnu-apps/NextGen/nextgensupport/hg18
elif [ $GENOME == "mm9" ]; then
  CHRSIZE=/nfs_exports/apps/gnu-apps/NextGen/nextgensupport/mm9/mm9.sizes
  ANNOT=/nfs_exports/apps/gnu-apps/NextGen/nextgensupport/mm9/chr2gene
else
  echo "wrong genome names, must be hg19, hg18 or mm9, case-sensitive"
  exit 1
fi

# Prepare commands file
SCRIPT=bw_presum.pl 
for st in "WG" "Exon" "CodingExon" "Intron" "Promoter"; do
  rm -f $SAMPLE.$st
  echo "#!/bin/bash" >$SAMPLE.$st
  echo "#BSUB -J $SAMPLE.$st" >>$SAMPLE.$st
  echo "#BSUB -e $SAMPLE.$st.err" >>$SAMPLE.$st
  echo "#BSUB -o $SAMPLE.$st.out" >>$SAMPLE.$st
  echo "hostname" >>$SAMPLE.$st
  echo "date" >>$SAMPLE.$st
  for i in `cat $CHRSIZE | cut -f1`; do
    if [ $i != "chrM" -a $i != "chrMT" ]; then
      CH=`echo $i | sed 's/chr//'`
      echo "$SCRIPT -chr $CH -cov_file $BW -annot_dir $ANNOT -chr_size $CHRSIZE -out_dir $COV_DIR/presum -sample $SAMPLE  -seg_type $st" \
      >>$SAMPLE.$st
    fi
  done
  echo "date" >>$SAMPLE.$st
done
exit;
