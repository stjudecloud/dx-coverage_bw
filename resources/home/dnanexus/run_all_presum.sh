#!/bin/bash
# Runs make_script_presum.sh to get a file with a list of commands and 
# runs the commands in that file.
#
# $1 = genome version, case sensitive (e.g. hg19, hg18, mm9)
# $2 = coverage directory, e.g.
#      /nfs_exports/genomes/1/PCGP/.allTumor/SJMB/NextGen/WholeGenome/Coverage
# $3 = subject ID, e.g., SJMB040
# $4 = sample name (pid only, not full barcode), e.g., SJMB040_D or SJMB040_G
# $5 = Segment types, e.g., "WG CodingExon Exon" for WGS

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

cmdsfile=`mktemp`
make_script_presum.sh "$@" > $cmdsfile
sh $cmdsfile
rm $cmdsfile 