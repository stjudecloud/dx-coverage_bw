#!/bin/bash
#Usage:  sh $COVSRC/auto_diff_track.sh G D SJxxx OUT_DIR CHRSIZE
#Assuming the denominator is SJxxx_G and numerator is SJxxx_D/R/X/... 
#Genomic average coverage is precomputed and saved in
#/nfs_exports/genomes/1/PCGP/BucketIntermediate/coverageWork/sampleinfo.txt
#SINFO=/nfs_exports/genomes/1/PCGP/BucketIntermediate/coverageWork/sampleinfo.txt
echo "running coverage diff..."
SA=$1 #GERMLINE SAMPLE
SB=$2 #TUMOR SAMPLE
#prefix of the output files ${SAMPLE}-${B}n.bw and ${SAMPLE}-diff_${B}n-${A}.bw
SAMPLE_PAIR=$3
OUT_DIR=$4 #this is /nfs_exports/genomes/1/projects/WHOLEGENOME/ClinicalPilot/.BucketFinal/SJEWS/Coverage/SJ001303 on research on clinical this is /clingen/dev/tartan/index/data/ClinicalPilot/ClinicalPilot/WHOLE_GENOME/coverage-paired

SINFOA=$5 #/nfs_exports/genomes/1/projects/WHOLEGENOME/ClinicalPilot/.BucketFinal/SJEWS/Coverage/Summary/SJEWS001303_D1.wg.cumulative.txt
SINFOB=$6 #/nfs_exports/genomes/1/projects/WHOLEGENOME/ClinicalPilot/.BucketFinal/SJEWS/Coverage/Summary/SJEWS001303_G1.wg.cumulative.txt
GENOME=$7
echo ${GENOME}


read case control < <(sn_string_to_pair.sh $SAMPLE_PAIR)
read case_subject case_disease case_type case_index case_sid < <(sn_parse.sh $case)
read control_subject control_disease control_type control_index control_sid < <(sn_parse.sh $control)

control_sample=`sn_build.sh $control_subject $control_disease $control_type $control_index`
case_sample=`sn_build.sh $case_subject $case_disease $case_type $case_index`

SAMPLE=${case_subject}


. import_config.sh genome $GENOME CHR_SIZES

if [ ! -f $OUT_DIR/$SA.bw ];then
echo "$OUT_DIR/$SA.bw does not exist!"
fi
if [ ! -f $OUT_DIR/$SB.bw ];then
echo "$OUT_DIR/$SB.bw does not exist!"
fi

#COVA=`cat $SINFO|perl -ne 's/\t/_/;print'|grep WholeGenome|grep $SA|cut -f7`
COVA=`cat ${SINFOA} | tail -1 | cut -f6` #20x coverage of wholegenome
#COVB=`cat $SINFO|perl -ne 's/\t/_/;print'|grep WholeGenome|grep $SB|cut -f7`
COVB=`cat ${SINFOB} | tail -1 | cut -f6` #20x coverage of wholegenome

if [ -z $COVA ]; then
echo "No mean covarage for $SA!"
exit 1
fi

if [ -z $COVB ]; then
echo "No mean covarage for $SB!"
exit 1
fi

date
echo "$SA coverage: $COVA"
echo "$SB coverage: $COVB"

#TEMP_DIR=`mktemp -d`
TEMP_DIR=`pwd`
cd $TEMP_DIR

echo "working directory: $TEMP_DIR"

while read i; do
	echo ${i};
	CHR=`echo -e "${i}" | cut -f1`
	SIZE=`echo -e "${i}" | cut -f2`

	echo "processing commands for $CHR:0-$SIZE"
	
	echo "bigWigToWig -chrom=$CHR -start=0 -end=${SIZE} $OUT_DIR/${SA}.bw ${TEMP_DIR}/${CHR}-A.wig" > ${CHR}.cmds
	echo "bigWigToWig -chrom=$CHR -start=0 -end=${SIZE} $OUT_DIR/${SB}.bw ${TEMP_DIR}/${CHR}-B.wig" >> ${CHR}.cmds
  
	#echo " normalize and substract..."
	echo "diff_track.pl ../${CHR}-A.wig ../${CHR}-B.wig $COVA $COVB $CHR"  >> ${CHR}.cmds
	echo "mv ${CHR}-diff.wig ../" >> ${CHR}.cmds
	echo "mv ${CHR}-Dn.wig ../" >> ${CHR}.cmds
	#echo " concatenate..";
done < $CHR_SIZES

ls chr*.cmds > diff.cmds

#submit job and wait
bsub_array_for_cmdfile.sh diff.cmds --log-dir `pwd`/logs -M 40000 -v 40000 -q clingen -J COV

while read i; do
	CHR=`echo -e "${i}" | cut -f1`
	COUNT=1
	NUMBER=0
	until [ "${NUMBER}" -eq "${COUNT}" ]
	do
		echo "Diff still running for ${CHR} $NUMBER / $COUNT"
		sleep 2
		NUMBER=`ls ${CHR}-diff.wig |wc -l`
	done
done < $CHR_SIZES

echo "diff tracks created for each chr, cat-ing and creating BW";

while read i; do

cat chr${i}-Dn.wig >> ${case}_normalized.wig
cat chr${i}-diff.wig >> ${SAMPLE_PAIR}_normalized-diff.wig

done < $CHR_SIZES

echo "wigToBigWig for ${B}n"
echo "wigToBigWig ${case}_normalized.wig $CHR_SIZES ${case_sample}_normalized.bw" > convert.cmds

echo "wigToBigWig for difference track" 
echo "wigToBigWig ${SAMPLE_PAIR}_normalized-diff.wig $CHR_SIZES ${SAMPLE_PAIR}_normalized-diff.bw" >> convert.cmds

bsub_array_for_cmdfile.sh convert.cmds --log-dir `pwd`/logs -M 40000 -v 40000 -q clingen -J convert

CHR=`echo -e "${i}" | cut -f1`

COUNT=4 #bw in same dir?
NUMBER=0
until [ "${NUMBER}" -eq "${COUNT}" ]
do
	echo "wigToBigWig still running for ${CHR} $NUMBER / $COUNT"
	sleep 5
	NUMBER=`ls *.bw |wc -l`
done

echo "copy files to $OUT_DIR"
#mv *.bw $OUT_DIR
#mv *.wig $OUT_DIR

#cd ../
#rm -fr $TEMP_DIR
echo "Done"
date
exit

