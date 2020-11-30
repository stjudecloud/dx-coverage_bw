#!/bin/bash
# Bambino 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://wiki.dnanexus.com/Developer-Portal for tutorials on how
# to modify this file.

main() {


  echo ""
  echo "=== Setup ==="
  echo "  [*] Downloading input files..." 
  echo "Value of BAM: '$BAM_path'"
  echo "Value of BAM index : '$BAM_INDEX_path'"
  dx-download-all-inputs --parallel > /dev/null

  ln -s $BAM_path $BAM_prefix.bam
  ln -s $BAM_INDEX_path $BAM_prefix.bam.bai


  ################
  # Housekeeping #
  ################

  echo "  [*] Performing some housekeeping..."
  # Fill in your application code here.
  bambino_path=`pwd`
  # Export JAVA classes
  export CLASSPATH=:$bambino_path/bambino-1.0.jar:$bambino_path/mysql-connector-java-5.1.38-bin.jar:$bambino_path/sam-1.65.jar:third_party.jar
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
  # Hg19 chromosome sizes

  echo "WARNING: Currently we are copying *.so.* files from a RHEL OS to /usr/lib64 so that sjWigToBigWig doesn't crash. This may be ill advised!!!"

  DEST_DIR=out
  CHR_SIZES="chrom.sizes.txt"
  dx download -o $CHR_SIZES project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/Homo_sapiens/${ref_name}/SUPPORT/canonical_chr_sizes.txt
  mkdir -p ${DEST_DIR}

  echo "  [*] Running Coverage ..." 
  # Get coverage using Ace2/Bambino SAMFinneyCoverage
  java -Xmx8192m Ace2.SAMFinneyCoverage -bam $BAM_prefix.bam -wig -stdout 1> ${DEST_DIR}/$BAM_prefix.wig

  echo "  [*] Making bw file ..." 
  #ldd $(which sjWigToBigWig)
  sjWigToBigWig ${DEST_DIR}/$BAM_prefix.wig $CHR_SIZES $BAM_prefix.bw

  echo "  [*] Finishing ..." 
  coverage_bw=$(dx upload $BAM_prefix.bw --brief)
  dx-jobutil-add-output coverage_bw "$coverage_bw" --class=file

}
