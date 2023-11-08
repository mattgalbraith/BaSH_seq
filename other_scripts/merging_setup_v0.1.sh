#!/bin/bash

SCRIPT_TITLE="merge_setup.sh"
SCRIPT_VERSION=0.1
# DATE: 11-29-2021
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to generate submit scriot for fastq merging. Requires tab-separated fastq_merging.txt file with 1) Desired SAMPLE_NAMES 2) read1/read2 3) full paths to fastq files to be merged in remaining columns.

# Run this script from newly created raw_merged_DD_MM_YYY dir

echo "Writing commands to $(pwd)/submit_fastq_merging.sh"

# write parameters to submit_fastq_merging.sh
echo "#!/bin/sh

SCRIPT_TITLE="$SCRIPT_TITLE"
SCRIPT_VERSION="$SCRIPT_VERSION"
# SUMMARY:
# This is an auto-generated script that contains commands to run merge fastq files specified in fastq_merging.txt.
# 
# This script should be run from the appropriate Project/raw_merged_MM_DD_YYYY directory as follows:
# sbatch submit_fastq_merging.sh
###### ^^ COPY AND PASTE COMMAND ABOVE TO RUN ^^ #######
########################################################
" > $(pwd)/submit_fastq_merging.sh

# write commands to submit_fastq_merging.sh
cat fastq_merging.txt | while IFS=$'\t' read F1 F2 F3 F4 F5
do 
	SUBMIT_LOG="$(pwd)/"$F1"_"$F2".log"
	echo "(nohup sh $HOME/PipeLinesScripts/other_scripts/FASTQ_MERGE_2.sh "$F1" "$F2" "$F3" "$F4" "$F5" &> "$SUBMIT_LOG") &"
done >> $(pwd)/submit_fastq_merging.sh

# need to execute commands above as subshells: (CMD) &
# then use "wait" which waits for child processes to complete
# then gather all std outs and errs to one log file:
echo "
wait

# Intialize log file:
echo  > $(pwd)/submit_fastq_merging.log

# Collecting individual std out and err files into one log file:" >> $(pwd)/submit_fastq_merging.sh

cat fastq_merging.txt | while IFS=$'\t' read F1 F2 F3 F4 F5
do
	SUBMIT_LOG="$(pwd)/"$F1"_"$F2".log"
	echo "cat "$SUBMIT_LOG" >> $(pwd)/submit_fastq_merging.log"
	echo "rm "$SUBMIT_LOG""
done >> $(pwd)/submit_fastq_merging.sh

echo "
# Intialize  file:
echo  > $(pwd)/sample_locations.txt

# Collecting individual sample locations into one  file:" >> $(pwd)/submit_fastq_merging.sh

cat fastq_merging.txt | while IFS=$'\t' read F1 F2 F3 F4 F5
do
	NAME="$(basename "$F3")"
	LOCATIONS_FILE=${NAME%.fastq.gz}_sample_locations.txt
	echo "cat "$LOCATIONS_FILE" >> $(pwd)/sample_locations.txt"
	echo "rm "$LOCATIONS_FILE""
done >> $(pwd)/submit_fastq_merging.sh

echo "Done."

exit