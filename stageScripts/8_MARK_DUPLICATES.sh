#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=8_MARK_DUPLICATES.sh
SCRIPT_VERSION=0.4
# DATE: 07-01-2021
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample (incl. both reads for PE) and is executed by <TTseq>_pipline.sh 
#
# v0.3: correcting $PICARD and command calls for Eureka
# Version 0.4 110923: updating to use PICARD Singularity/Apptainer on Proton2

# variables from command line via <TTseq>_pipline.sh:
SAMPLE_DIR=${1}
ALIGNMENT_DIRNAME=${2}
BAM_IN_FILENAME=${3}
BAM_OUT_FILENAME=${4}
PICARD_MEM=${5}
DUPLICATES=${6}
PICARD_SIF=${7}
THIS_ANALYSIS_DIR=${8}
#
# PICARD=$(realpath $(which picard)) # work around for conda install not working with java -jar call
#
# other variables:
# BAM_IN_FILE (see below)
# BAM_OUT_FILE (see below)
# TMP_DIR (see below; $TMPDIR is built in to Rosalind profile)
# MARKDUPS_VERSION="$(java -Xmx1G -jar $PICARD.jar MarkDuplicates --version 2>&1)"
MARKDUPS_VERSION=`singularity run "$PICARD_SIF" java -jar /picard.jar MarkDuplicates --version 2>&1`

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for "$SCRIPT_TITLE":
(1) SAMPLE_DIR: "$SAMPLE_DIR"
(2) ALIGNMENT_DIRNAME: "$ALIGNMENT_DIRNAME"
(3) BAM_IN_FILENAME: "$BAM_IN_FILENAME"
(4) BAM_OUT_FILENAME: "$BAM_OUT_FILENAME"
(5) PICARD_MEM: "$PICARD_MEM"
(6) DUPLICATES: "$DUPLICATES"
(7) PICARD_SIF: "$PICARD_SIF"
(8) THIS_ANALYSIS_DIR: "$THIS_ANALYSIS_DIR"
Picard MarkDuplicates version: "$MARKDUPS_VERSION"
PICARD MarkDuplicates options:
ASSUME_SORTED=true
CREATE_INDEX=true
TMP_DIR=$TMPDIR
MAX_RECORDS_IN_RAM=1000000
"

EXPECTED_ARGS=8
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo "Usage: "$SCRIPT_TITLE" <SAMPLE_DIR> <ALIGNMENT_DIRNAME> <BAM_IN_FILENAME> <BAM_OUT_FILENAME> <PICARD_MEM> <mark/remove duplicates> <PICARD_SIF> <THIS_ANALYSIS_DIR>
        ${red}ERROR - expecting "$EXPECTED_ARGS" ARGS but "$#" were provided:${NC}
        "$@"
        "
        exit 1
fi

# Check for BAM file in alignment directory
BAM_IN_FILE=$SAMPLE_DIR/$ALIGNMENT_DIRNAME/$BAM_IN_FILENAME
if [ ! -f $BAM_IN_FILE ]
then
	echo "${red}ERROR - file does not exist: "$BAM_IN_FILE"
	"
	exit 1
fi

if [ "$DUPLICATES" = "remove" ]
then
	echo "removing duplicates
	"
	REMOVE_DUPLICATES="true"
else
	echo "marking duplicates
	"
	REMOVE_DUPLICATES="false"
fi

BAM_OUT_FILE="$SAMPLE_DIR"/"$ALIGNMENT_DIRNAME"/"$BAM_OUT_FILENAME"

echo -e "
${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"

	singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$PICARD_SIF" java -XX:ParallelGCThreads=5 -Xmx"$PICARD_MEM" -jar /picard.jar MarkDuplicates \
		INPUT="$BAM_IN_FILE" \
		OUTPUT="$BAM_OUT_FILE" \
		METRICS_FILE="$BAM_IN_FILE".duplicate_metrics.txt \
		REMOVE_DUPLICATES="$REMOVE_DUPLICATES" \
		ASSUME_SORTED=true \
		CREATE_INDEX=true \
		TMP_DIR=$TMPDIR/"$SLURM_JOB_ID" \
		MAX_RECORDS_IN_RAM=1000000

	# check output status
	if [ $? -ne 0 ]
	then
		echo -e "${red}PICARD MarkDuplicates failed
		"
		exit 1
	fi

	# now copy the .bai file to .bam.bai file (some tools need it in that format)
	BAI=`echo $BAM_OUT_FILE | sed 's/\.bam/\.bai/'`
	cp $BAI $BAM_OUT_FILE.bai

echo -e "${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]${NC}"

