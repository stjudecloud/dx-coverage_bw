#!/usr/bin/env bash
#
# Moves bigwig files to the given destination by first copying it
# and then removing the source. It also creates a symlink from the source
# to the destination.
#
# This script is useful to move the bigwig files from their original
# location to a location where the genome browser can access them.
#
# Params:
# $1 = Path to bigwig file
# $2 = Path to the destination directory where the bigwig file needs to be moved to
#

# Show usage information if no parameters were sent
if [ "$#" -lt 2 ]; then about.sh $0; exit 1; fi

# Get params
bigwig=$1
destination=$2

# If bigwig is already a symlink, then we do not need to proceed
if [ -h "$bigwig" ]
then
  echo "Note: Bigwig file $bigwig is already a symlink, so no need to move and symlink."
  exit
fi

# Get bigwig filename
bigwig_filename=`basename $bigwig`

# Copy bigwig to destination
cp $bigwig ${destination}/
exit_code=$?
if [ "$exit_code" -ne 0 ]
then
  echo "Error: There was an error copying $bigwig to $destination. Exiting..."
  exit $exit_code
fi

# Remove original bigwig
rm $bigwig
exit_code=$?
if [ "$exit_code" -ne 0 ]
then
  echo "Error: There was an error removing bigwig file $bigwig. Exiting..."
  exit $exit_code
fi

# Create symlink from original location to new location
ln -s ${destination}/${bigwig_filename} $bigwig
exit_code=$?
if [ "$exit_code" -ne 0 ]
then
  echo "Error: There was an error creating a symbolic link from $bigwig to ${destination}/${bigwig_filename}. Exiting..."
  exit $exit_code
fi
