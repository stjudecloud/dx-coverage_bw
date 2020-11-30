#/usr/bin/sh


run_dir=`tartan_run.pl new type=GOLD_GENE_COVERAGE_CHECK`
echo "TARTAN DIR: ${run_dir}"

#PAIR=$1 #SJLGG022_D_SJLGG022_G
case_sample=$1
control_sample=$2
CUTOFF=$3 #70
ROOT=$4 #/clingen/dev/tartan/index/data/ClinicalPilot/ClinicalPilot/
GENOME=$5

#new location for big wig files
#/clingen/dev/tartan/index/data/ClinicalPilot/ClinicalPilot/${case_sample}/WHOLE_GENOME/coverage/${case_sample}.bw


echo ${GENOME}

#chr2gene=$7 #/nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/chr2gene

. import_config.sh genome $GENOME CHR2GENE CLINCLS_GOLD_GENE_LIST_FILE CLINCLS_GL_GOLD_NON_CANCER_GENE_FILE CLINCLS_GL_PCGP_GOLD_FILE GENE_TRANSCRIPT_MATRIX CLINCLS_CANCER_RELATED_GENES_FILE

#read case control < <(sn_string_to_pair.sh $PAIR)
#read case_subject case_disease case_type case_index case_sid < <(sn_parse.sh $case)
#read control_subject control_disease control_type control_index control_sid < <(sn_parse.sh $control)

#control_sample=`sn_build.sh $control_subject $control_disease $control_type $control_index`
#case_sample=`sn_build.sh $case_subject $case_disease $case_type $case_index`	

#WGS_D_BW=${BW_WGS_ROOT}/${case_sample}.bw
WGS_D_BW=`tartan_run.pl $run_dir addInput src=${ROOT}/${case_sample}/WHOLE_GENOME/coverage/${case_sample}.bw dest=WGS/${case_sample}.bw`
#WGS_G_BW=${BW_WGS_ROOT}/${control_sample}.bw
WGS_G_BW=`tartan_run.pl $run_dir addInput src=${ROOT}/${control_sample}/WHOLE_GENOME/coverage/${control_sample}.bw dest=WGS/${control_sample}.bw`
#EXCAP_D_BW=${BW_EXCAP_ROOT}/${case_sample}.bw
EXCAP_D_BW=`tartan_run.pl $run_dir addInput src=${ROOT}/${case_sample}/EXOME/coverage/${case_sample}.bw dest=EXCAP/${case_sample}.bw`
#EXCAP_G_BW=${BW_EXCAP_ROOT}/${control_sample}.bw
EXCAP_G_BW=`tartan_run.pl $run_dir addInput src=${ROOT}/${control_sample}/EXOME/coverage/${control_sample}.bw dest=EXCAP/${control_sample}.bw`

tartan_run.pl $run_dir startWork

cd `tartan_run.pl $run_dir workspaceDir`


#SOMATIC_GOLD=`grep -w "\-gold" ${MEDAL_CONFIG} | cut -f2 -d" "`
SOMATIC_GOLD=${CLINCLS_GOLD_GENE_LIST_FILE}
#GERMLINE_GOLD_CANCER=`grep -w "\-gl-gold-cancer" ${MEDAL_CONFIG} | cut -f2 -d" "`
cat ${CLINCLS_CANCER_RELATED_GENES_FILE} | cut -f1 > GERMLINE_GOLD_CANCER.tmp
GERMLINE_GOLD_CANCER=GERMLINE_GOLD_CANCER.tmp
#GERMLINE_GOLD_PCGP=`grep -w "\-gl-pcgp-gold-db"  ${MEDAL_CONFIG} | cut -f2 -d" "`
GERMLINE_GOLD_PCGP=${CLINCLS_GL_PCGP_GOLD_FILE}
#GERMLINE_GOLD_NON_CANCER=`grep -w "\-gl-gold-non-cancer" ${MEDAL_CONFIG} | cut -f2 -d" "`
GERMLINE_GOLD_NON_CANCER=${CLINCLS_GL_GOLD_NON_CANCER_GENE_FILE}

ISOFORM_LIST=${GENE_TRANSCRIPT_MATRIX}

#MEDAL_CONFIG=`tartan_run.pl $run_dir addInput src=${MEDAL_CONFIG} dest=.`
#SOMATIC_GOLD=`tartan_run.pl $run_dir addInput src=${SOMATIC_GOLD} dest=.`
#GERMLINE_GOLD_CANCER=`tartan_run.pl $run_dir addInput src=${GERMLINE_GOLD_CANCER} dest=.`
#GERMLINE_GOLD_PCGP=`tartan_run.pl $run_dir addInput src=${GERMLINE_GOLD_PCGP} dest=.`

#ISOFORM_LIST=`tartan_run.pl $run_dir addInput src=${ISOFORM_LIST} dest=.`

echo -e "${WGS_D_BW}\n${EXCAP_D_BW}" > ${case_sample}_bw.lst
echo -e "${WGS_G_BW}\n${EXCAP_G_BW}" > ${control_sample}_bw.lst

#echo "SOMATIC_GOLD = '${SOMATIC_GOLD}'
#GERMLINE_GOLD_CANCER = '${GERMLINE_GOLD_CANCER}'
#GERMLINE_GOLD_PCGP = '${GERMLINE_GOLD_PCGP}'"

GERMLINE_GOLD_PCGP_TMP=GERMLINE_GOLD_PCGP.tmp
cut -f1 ${GERMLINE_GOLD_PCGP} | grep -v GeneName | sort | uniq > ${GERMLINE_GOLD_PCGP_TMP}


GENES[0]=`echo "SOMATIC_GOLD;${SOMATIC_GOLD}"`
GENES[1]=`echo "GERMLINE_GOLD_CANCER;${GERMLINE_GOLD_CANCER}"`
GENES[2]=`echo "GERMLINE_GOLD_PCGP;${GERMLINE_GOLD_PCGP_TMP}"`
GENES[3]=`echo "GERMLINE_GOLD_NON_CANCER;${GERMLINE_GOLD_NON_CANCER}"`


for line in "${GENES[@]}"; do
	echo ${line}
	GENE_LIST_NAME=`echo ${line} | cut -f1 -d";"`
	GENE_LIST=`echo ${line} | cut -f2 -d";"`
	ISOFORM_LIST_OUT=${GENE_LIST_NAME}.isoforms
	for i in `cat ${GENE_LIST}`; do
		GENE_NAME=`echo -e "${i}" | cut -f1`
		if [ "${GENE_NAME}" != "GeneName" ]
		then
			#echo -e "grep -w ${GENE_NAME} ${ISOFORM_LIST} | head -1 | cut -f3 | sed 's/^M//g' >> ${ISOFORM_LIST_OUT}"
			grep -w ${GENE_NAME} ${ISOFORM_LIST} | head -1 | cut -f3 | sed 's/^M//g' >> ${ISOFORM_LIST_OUT}
		fi
	done
done

echo "gene list created"

FINAL_ISOFORMS='isoforms_to_check.txt'

echo "isoform list created"

cat *.isoforms | sort | uniq > ${FINAL_ISOFORMS}
echo "submitting jobs"

bsub -q clingen -M 8000 -v 8000 -J ${control_sample} -o ${control_sample}.cluster.out -e ${control_sample}.cluster.err bw_exon_coverage.pl -isoform-list ${FINAL_ISOFORMS} -chr2gene ${CHR2GENE} -bw-list ${control_sample}_bw.lst -suffix combined -combine

bsub -q clingen -M 8000 -v 8000 -J ${case_sample} -o ${case_sample}.cluster.out -e ${case_sample}.cluster.err bw_exon_coverage.pl -isoform-list ${FINAL_ISOFORMS} -chr2gene ${CHR2GENE} -bw-list ${case_sample}_bw.lst -suffix combined -combine


#wait unitl jobs finnish

NUMBER=0
COUNT=2
until [ "${NUMBER}" -eq "${COUNT}" ]
do
		echo "Coverage Still Running $NUMBER / $COUNT"
		sleep 5
		NUMBER=`ls *.combined.coverage.tab | wc -l`
done



#write config file for perl script
CONFIG_D=`pwd`/${case_sample}_gene_coverage_check.config
CONFIG_G=`pwd`/${control_sample}_gene_coverage_check.config
OUTPUT=`pwd`/${case_sample}_gene_coverage_report.tab
PWD=`pwd`
echo "SAMPLE	${case_sample}
COVERAGE_FILE	${PWD}/${case_sample}.bw.combined.coverage.tab
OUTPUT_FILE	${OUTPUT}
CUTOFF	${CUTOFF}" > $CONFIG_D

echo "SAMPLE	${control_sample}
COVERAGE_FILE	${control_sample}.bw.combined.coverage.tab
OUTPUT_FILE	${control_sample}_gene_coverage_report.tab
CUTOFF	${CUTOFF}" > $CONFIG_G

#now run perl script
echo "RUNNING: coverage_summary_gene.pl -c $CONFIG_D"
coverage_summary_gene.pl -c $CONFIG_D
echo "RUNNING: coverage_summary_gene.pl -c $CONFIG_G"
coverage_summary_gene.pl -c $CONFIG_G

mv *coverage.tab ../output/
mv *coverage_report.tab ../output/

tartan_run.pl $run_dir endWork
