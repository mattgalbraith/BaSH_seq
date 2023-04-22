#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=5_ADD_RGID.sh
SCRIPT_VERSION=0.3
# DATE: 07-01-2021
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample (incl. both reads for PE) and is executed by <TTseq>_pipline.sh
#
# v0.3: correcting $PICARD and command calls for Eureka

# variables from command line via <TTseq>_pipline.sh:
SAMPLE_NAME=$1
PLATFORM=$2
DATE=$3
PI=$4
LIB=$5
SEQ_CORE=$6
SEQ_ID=$7
EXPERIMENT=$8
SAMPLE_DIR=${9}
ALIGNMENT_DIRNAME=${10}
BAM_IN_FILENAME=${11}
BAM_OUT_FILENAME=${12}
PICARD_MEM=${13}
#
PICARD=$(realpath $(which picard)) # work around for conda install not working with java -jar call
#
ADD_RGID_VERSION="$(java -Xmx1G -jar $PICARD.jar AddOrReplaceReadGroups --version 2>&1)"

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for "$SCRIPT_TITLE":
(1) SAMPLE_NAME: "$SAMPLE_NAME"
(2) PLATFORM: "$PLATFORM"
(3) DATE: "$DATE"
(4) PI: "$PI"
(5) LIBRARY: "$LIB"
(6) SEQ_CORE: "$SEQ_CORE"
(7) SEQ_ID: "$SEQ_ID"
(8) EXPERIMENT: "$EXPERIMENT"
(9) SAMPLE_DIR: "$SAMPLE_DIR"
(10) ALIGNMENT_DIRNAME: "$ALIGNMENT_DIRNAME"
(11) BAM_IN_FILENAME: "$BAM_IN_FILENAME"
(12) BAM_OUT_FILENAME: "$BAM_OUT_FILENAME"
(13) PICARD_MEM: "$PICARD_MEM"
Picard AddOrReplaceReadGroups version: "$ADD_RGID_VERSION"
Picard AddOrReplaceReadGroups options:
CREATE_INDEX=true
TMP_DIR="$TMPDIR"
MAX_RECORDS_IN_RAM=1000000
"

EXPECTED_ARGS=13
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo "Usage: "$SCRIPT_TITLE" <SAMPLE_NAME> <PLATFORM> <DATE> <PI (investigator code)> <LIBRARY (PE/SE)> <SEQ_CORE> <SEQ_ID> <EXPERIMENT> <SAMPLE_DIR> <ALIGNMENT_DIRNAME> <BAM_IN_FILENAME> <BAM_OUT_FILENAME> <PICARD_MEM>
        ${red}ERROR - expecting "$EXPECTED_ARGS" ARGS but "$#" were provided:${NC}
        "$@"
        "
        exit 1
fi

# Check for BAM file in alignment directory
if [ ! -f $SAMPLE_DIR/$ALIGNMENT_DIRNAME/$BAM_IN_FILENAME ]
then
	echo "${red}ERROR - file does not exist: $SAMPLE_DIR/$ALIGNMENT_DIRNAME/$BAM_IN_FILENAME
	"
	exit 1
fi

echo -e "
${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]${NC}"

	srun java -Xmx"$PICARD_MEM" -jar $PICARD.jar AddOrReplaceReadGroups \
		INPUT="$SAMPLE_DIR"/"$ALIGNMENT_DIRNAME"/"$BAM_IN_FILENAME" \
		OUTPUT="$SAMPLE_DIR"/"$ALIGNMENT_DIRNAME"/"$BAM_OUT_FILENAME" \
		TMP_DIR=$TMPDIR \
		RGID="$SAMPLE_NAME" \
		RGLB="$SAMPLE_NAME" \
		RGPL="$PLATFORM" \
		RGPU="$SAMPLE_NAME" \
		RGSM="$SAMPLE_NAME" \
		RGCN="$SEQ_CORE" \
		RGDT="$DATE" \
		RGDS="$EXPERIMENT"-"$PI"-"$LIB" \
		CREATE_INDEX=true \
		MAX_RECORDS_IN_RAM=1000000

	# check output status
	if [ $? -ne 0 ]
	then
			echo "PICARD AddOrReplaceReadGroups failed"
			exit 1
	fi

echo -e "${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]${NC}"

