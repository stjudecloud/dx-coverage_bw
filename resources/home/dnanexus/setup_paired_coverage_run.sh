#!/usr/bin/env bash
# Writes scripts for paired-coverage runs.
#
# $1 = genome
# $2 = config file
# $3 = run directory
# $4 = data (input) directory
# $5 = output directory

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
GENOME=$1
ANLS_CONFIG=$2
RUN_DIR=$3
DATA_DIR=$4
OUTPUT_DIR=$5

# Get the anls config if not specified
if [ "$ANLS_CONFIG" == "" ]; then ANLS_CONFIG=`ls $RUN_DIR/config*`; fi

# Validate anls config
if [ ! -f "$ANLS_CONFIG" ]
then echo "Analysis config file not found: $ANLS_CONFIG" >&2; exit 1
fi

# Load configs
echo
echo "Reading config files (you may see some error messages--these are OK):"
echo "* Genome level"
echo "* App level"
echo "..."
. import_config.sh genome $GENOME
. import_config.sh app coverage-paired

# Load step library
. steplib.sh

# Create run directory if it does not exist
if [ ! -d "$RUN_DIR" ]; then mkdir -p "$RUN_DIR"; fi

# Copy in analysis file
cp -v $ANLS_CONFIG $RUN_DIR/

# Set the run directory
cd $RUN_DIR
set_step_script_dir $RUN_DIR

echo "Setup paired-coverage run"

#######################################
# Step: Paired-coverage
#######################################
# Initialize step
init_step cov-paired

# Make-commands script
cat > `get_step_make_cmds_script` <<EOF
#!/usr/bin/env bash
# Write paired-coverage commands to commands file
cat $ANLS_CONFIG | while read bam_pair case_bam control_bam
do 
  for gender in M F
  do
    echo covsum_paired_wg.sh \$case_bam \$control_bam \$gender $DATA_DIR/\$bam_pair
  done
done > $RUN_DIR/cmds-01.sh
EOF

# Submit script
write_step_submit_script

#######################################
# Step: Bigwig to wig
#######################################
# Initialize step
init_step cov-bw2wig

# Make-commands script
cat > `get_step_make_cmds_script` <<EOF
#!/usr/bin/env bash
# Write bigwig-to-wig commands to commands file
cat $ANLS_CONFIG | while read bam_pair case_bam control_bam
do
  mkdir -p $DATA_DIR/\$bam_pair/wig_files
  make_script_bigwig2wig_per_chr.sh $GENOME \$case_bam $DATA_DIR/\$bam_pair/Summary $DATA_DIR/\$bam_pair/wig_files
  exit_code=\$?
  if [ "\$exit_code" -ne 0 ]; then exit \$exit_code; fi
  make_script_bigwig2wig_per_chr.sh $GENOME \$control_bam $DATA_DIR/\$bam_pair/Summary $DATA_DIR/\$bam_pair/wig_files
  exit_code=\$?
  if [ "\$exit_code" -ne 0 ]; then exit \$exit_code; fi
done > $RUN_DIR/cmds-02.sh
EOF

# Submit script
write_step_submit_script

#######################################
# Step: Diff
#######################################
# Initialize step
init_step cov-diff

# Make-commands script
cat > `get_step_make_cmds_script` <<EOF
#!/usr/bin/env bash
# Write coverage-diff commands to commands file
cat $ANLS_CONFIG | while read bam_pair case_bam control_bam
do 
  make_script_diff_track.sh $GENOME \$case_bam \$control_bam $DATA_DIR/\$bam_pair/wig_files $DATA_DIR/\$bam_pair/wig_files $DATA_DIR/\$bam_pair/Summary/\${case_bam}.coverage.summary.table.txt
  exit_code=\$?
  if [ "\$exit_code" -ne 0 ]; then exit \$exit_code; fi
done > $RUN_DIR/cmds-03.sh
EOF

# Submit script
write_step_submit_script

#######################################
# Step: Normalized-diff
#######################################
# Initialize step
init_step cov-normalized-diff

# Make-commands script
cat > `get_step_make_cmds_script` <<EOF
#!/usr/bin/env bash
# Write normalized-coverage-diff commands to commands file

# Write normalized-coverage-diff commands to commands file
cat $ANLS_CONFIG | while read bam_pair case_bam control_bam
do
  # Set input and output dir
  input_dir=$DATA_DIR/\$bam_pair/wig_files
  output_dir=$DATA_DIR/\$bam_pair/Summary

  # Concatenate per-chr wig results
  cat \$input_dir/chr*-Dn.wig > \$output_dir/\${case_bam}.normalized.wig
  cat \$input_dir/chr*-diff.wig > \$output_dir/\${bam_pair}.normalized-diff.wig

  # Create command to convert from wig to bigwig
  echo "wigToBigWig \$output_dir/\${case_bam}.normalized.wig $CHR_SIZES \$output_dir/\${case_bam}.normalized.bw" 
  echo "wigToBigWig \$output_dir/\${bam_pair}.normalized-diff.wig $CHR_SIZES \$output_dir/\${bam_pair}.normalized-diff.bw"
done > $RUN_DIR/cmds-04.sh
EOF

# Submit script
write_step_submit_script

#######################################
# Step: Move results to output directory
#######################################
# Initialize step
init_step cp-results

# Local-work script
cat > `get_step_local_work_script` <<EOF
#!/usr/bin/env bash
# Write script for copying results to output directory
while read bam_pair case_bam control_bam
do 
  # Create output directory if it doesn't exist
  mkdir -p $OUTPUT_DIR/\$bam_pair

  # Move important result files to output directory
  result_dir=$DATA_DIR/\$bam_pair/Summary
  for i in \$result_dir/*.png \$result_dir/*.coverage.summary.table.txt \$result_dir/*normalized.bw \$result_dir/*normalized-diff.bw
  do
    # Get file name from the path
    f=\`echo \$i | awk -F/ '{print \$NF}'\`

    # Move result file
    mv \$i $OUTPUT_DIR/\$bam_pair/\$f
  done

  # Move important files in the 'F' subdir in the input dir to the 'F' subdir in the output dir
  mkdir -p $OUTPUT_DIR/\$bam_pair/F
  for i in \$result_dir/F/*.png \$result_dir/F/*.coverage.summary.table.txt
  do
    # Get file name from the path
    f=\`echo \$i | awk -F/ '{print \$NF}'\`

    # Move result file
    mv \$i $OUTPUT_DIR/\$bam_pair/F/\$f
  done
done < $ANLS_CONFIG
EOF
