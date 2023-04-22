#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=7_SORT_BAM.sh
SCRIPT_VERSION=0.3
# DATE: 07-01-2021
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample (incl. both reads for PE) and is executed by <TTseq>_pipline.sh 
#
# v0.3: correcting $PICARD and command calls for Eureka

# variables from command line via <TTseq>_pipline.sh:
SAMPLE_DIR=$1
ALIGNMENT_DIRNAME=$2
BAM_IN_FILENAME=$3
BAM_OUT_FILENAME=$4
PICARD_MEM=$5
#
PICARD=$(realpath $(which picard)) # work around for conda install not working with java -jar call
#
# Other variables:
SORTSAM_VERSION="$(java -Xmx1G -jar $PICARD.jar SortSam --version 2>&1)"

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for SortBam.sh:
(1) SAMPLE_DIR: "$SAMPLE_DIR"
(2) ALIGNMENT_DIRNAME: "$ALIGNMENT_DIRNAME"
(3) BAM_IN_FILENAME: "$BAM_IN_FILENAME"
(4) BAM_OUT_FILENAME: "$BAM_OUT_FILENAME"
(5) PICARD_MEM: "$PICARD_MEM"
PICARD SortSam options:
SORT_ORDER=coordinate
CREATE_INDEX=true
TMP_DIR="$TMPDIR"
MAX_RECORDS_IN_RAM=1000000
Picard SortSam version: "$SORTSAM_VERSION"
Picard SortSam options:
SORT_ORDER=coordinate
CREATE_INDEX=true
MAX_RECORDS_IN_RAM=1000000
"

EXPECTED_ARGS=5
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo "Usage: "$SCRIPT_TITLE" <SAMPLE_DIR> <ALIGNMENT_DIRNAME> <BAM_IN_FILENAME> <BAM_OUT_FILENAME> <PICARD_MEM>
        ${red}ERROR - expecting "$EXPECTED_ARGS" ARGS but "$#" were provided:${NC}
        "$@"
        "
        exit 1
fi

# Check for BAM file in alignment directory
BAM_IN_FILE=$SAMPLE_DIR/$ALIGNMENT_DIRNAME/$BAM_IN_FILENAME 
if [ ! -f $BAM_IN_FILE ]
then
	echo "${red}ERROR - file does not exist: $BAM_IN_FILE"
	echo
	exit 1
fi

BAM_OUT_FILE=$SAMPLE_DIR/$ALIGNMENT_DIRNAME/$BAM_OUT_FILENAME

echo -e "${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}
"
	srun java -XX:ParallelGCThreads=5 -Xmx"$PICARD_MEM" -jar $PICARD.jar SortSam \
		INPUT="$BAM_IN_FILE" \
		OUTPUT="$BAM_OUT_FILE" \
		SORT_ORDER=coordinate \
		CREATE_INDEX=true \
		TMP_DIR="$TMPDIR" \
		MAX_RECORDS_IN_RAM=1000000

	# check output status
	if [ $? -ne 0 ]
	then
		echo "${red}PICARD SortSam failed"
		exit 1
	fi

echo -e "${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]${NC}"

