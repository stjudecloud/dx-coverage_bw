#!/bin/bash
# Writes step scripts for coverage runs.
#
# $1 = genome
# $2 = coverage configuration file 
# $3 = run directory
# $4 = input directory
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
echo "* Application level"
echo "..."
. import_config.sh genome $GENOME CHR_SIZES
. import_config.sh app coverage

# Read step script
. steplib.sh

# Create run directory if it does not exist
if [ ! -d "$RUN_DIR" ]; then mkdir -p "$RUN_DIR"; fi

# Copy in analysis file
cp -v $ANLS_CONFIG $RUN_DIR/

# Set the run directory
cd $RUN_DIR
set_step_script_dir $RUN_DIR

echo "Setting up coverage run..."

# Get the list of samples on a single line for looping
SAMPLES=`cat $ANLS_CONFIG | tr '\n' ' '`

#######################################
# Step: Create coverage bigwig
#######################################
# Initialize step
init_step cov-bw

# Create cmds-script
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
# This script writes coverage commands to the commands file

# Write coverage command for each sample
for sample in $SAMPLES
do 
  echo coverage_bw.sh $DATA_DIR/\$sample/\$sample.bam $CHR_SIZES $OUTPUT_DIR/\$sample.bw
done > `get_step_cmds_file`
EOF

# Create submit-script
write_step_submit_script

#######################################
# End of Steps
#######################################

# Print run dir
echo "...done"
echo "Run dir is $RUN_DIR"
