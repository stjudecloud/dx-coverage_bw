#!/bin/bash
# Writes scripts for coverage-post runs.
#
# $1 = target
# $2 = genome
# $3 = config file
# $4 = run directory
# $5 = data (input) directory
# $6 = output directory

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
TARGET=$1
GENOME=$2
ANLS_CONFIG=$3
RUN_DIR=$4
DATA_DIR=$5
OUTPUT_DIR=$6

# Get the anls config if not specified
if [ "$ANLS_CONFIG" == "" ]; then ANLS_CONFIG=`ls $RUN_DIR/config*`; fi

# Validate anls config
if [ ! -f "$ANLS_CONFIG" ]
then echo "Analysis config file not found: $ANLS_CONFIG" >&2; exit 1
fi

# Load configs
echo
echo "Reading config files (you may see some error messages--these are OK):"
echo "* Application level"
echo "* Sequencing target"
echo "..."
. import_config.sh app coverage-post
. import_config.sh target $TARGET SEGMENT_TYPES PRIMARY_SEGMENT_TYPE MIN_COVGT20X_CASE MIN_COVGT20X_CONTROL COV_THRESHOLD

if [ "$COV_THRESHOLD" == "" ] 
then 
  COV_THRESHOLD=20
fi 

# Load step library
. steplib.sh

# Create run directory if it does not exist
if [ ! -d "$RUN_DIR" ]; then mkdir -p "$RUN_DIR"; fi

# Copy in analysis file
cp -v $ANLS_CONFIG $RUN_DIR/

# Create some dirs for status under the run dir
# These are for small files that we want to keep but not index
outfile_case=$RUN_DIR"/coverage_qc_stats_case"
mkdir -p $outfile_case
outfile_control=$RUN_DIR"/coverage_qc_stats_control"
mkdir -p $outfile_control
cov_val_dir=$RUN_DIR"/coverage_cumulative_stats"
mkdir -p $cov_val_dir

# Set the run directory
cd $RUN_DIR
set_step_script_dir $RUN_DIR

echo "Setup coverage-post run"

# Get the list of samples on a single line for looping, and save to samples.txt
SAMPLES=`cut -f1 $ANLS_CONFIG | tee samples.txt | tr '\n' ' '`

#######################################
# Step: Presum
#######################################
# Initialize step
init_step cov-presum

# Make-commands script
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
# Write presum-compute commands to commands file
for sample in $SAMPLES
do make_script_presum.sh $GENOME $DATA_DIR/\$sample/coverage \$sample \$sample $SEGMENT_TYPES
done > `get_step_cmds_file`
EOF

# Submit script
write_step_submit_script

#######################################
# Step: Summary
#######################################
# Initialize step
init_step cov-summary

# Make-commands script
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
# Write coverage-summary commands to commands file
for sample in $SAMPLES
do
  for gender in M F
  do echo covsum.sh \$sample \$gender $DATA_DIR/\$sample/coverage $GENOME $SEGMENT_TYPES
  done
done > `get_step_cmds_file`
EOF

# Submit script
write_step_submit_script

#######################################
# Step: Local stats and plotting
#######################################
if [ "$PRIMARY_SEGMENT_TYPE" ]
then
  # Initialize step
  init_step cov-qc

  # Segment type is lowercase in the filename
  segtypepart=`echo $PRIMARY_SEGMENT_TYPE | tr '[:upper:]' '[:lower:]'`
  
  # Run QC stat gathering code locally, since it is fast.
  cat > `get_step_qc_script` <<EOF
#!/bin/bash
# Gather coverage data for QC
. import_config.sh target $TARGET MIN_COVGT${COV_THRESHOLD}X_CASE MIN_COVGT${COV_THRESHOLD}X_CONTROL
echo -e "# Sample\tStatus\tCoverageValue\tCutoff" > $OUTPUT_DIR/bam_qc_summary.txt
cat $ANLS_CONFIG | while read sample case_or_control
do
  # Get ${COV_THRESHOLD}x coverage value 
  file=$DATA_DIR/\$sample/coverage/Summary/\$sample.$segtypepart.cumulative.txt

  qc_coverage.pl -t $TARGET -f \$file -s \$sample -c \$case_or_control
done >> $OUTPUT_DIR/bam_qc_summary.txt

#Extract coverage values for plots
echo -e "Sample\t10x\t20x\t30x\t40x\t45x" > $OUTPUT_DIR/cumulative_coverage.txt
for sample in $SAMPLES
 do
  grep "rcumsum" $DATA_DIR/\$sample/coverage/Summary/\$sample.$segtypepart.cumulative.txt | cut -f5,6,7,8,9 | awk -v s=\$sample '{ print s, \$0 }'
done >> $OUTPUT_DIR/cumulative_coverage.txt

if [ "$TARGET" == "WHOLE_GENOME" ]
then 
for sample in $SAMPLES
 do
   cd $DATA_DIR/\$sample/coverage/Summary; unpaired_cov_hist.R \$sample 101
done
fi
EOF

  cat > `get_step_local_work_script` <<EOF
#!/bin/bash
# Move files from data directory to output directory
cat $ANLS_CONFIG | while read sample case_or_control
do
  mkdir -p $OUTPUT_DIR/\$sample/F
  for i in $DATA_DIR/\$sample/coverage/Summary/*.*
  do
    f=\`echo \$i | awk -F/ '{print \$NF}'\`
    mv \$i $OUTPUT_DIR/\$sample/\$f
  done
  
  for i in $DATA_DIR/\$sample/coverage/Summary/F/*.*
  do
    f=\`echo \$i | awk -F/ '{print \$NF}'\`
    mv \$i $OUTPUT_DIR/\$sample/F/\$f
  done
done
EOF
fi

#######################################
# End of Steps
#######################################

# Print instructions
echo "Run dir is $RUN_DIR"
