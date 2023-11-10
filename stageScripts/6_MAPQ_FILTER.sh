#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=6_MAPQ_FILTER.sh
SCRIPT_VERSION=0.3
# DATE: 06-21-2018
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample (incl. both reads for PE) and is executed by <XXseq>_pipeline.sh 
#
# Version 0.3 110923: updating to use PICARD Singularity/Apptainer on Proton2

# variables from command line via <XXseq>_pipeline.sh:
SAMPLE_DIR=${1}
ALIGNMENT_DIRNAME=${2}
BAM_IN_FILENAME=${3}
BAM_OUT_FILENAME=${4}
MIN_MAPQ=${5}
SAMTOOLS_SIF=${6}
THIS_ANALYSIS_DIR=${7}
# Other variables:
# SAMTOOLS_VERSION="$(samtools --version 2>&1)"
SAMTOOLS_VERSION=`singularity run "$SAMTOOLS_SIF" samtools --version | head -n1`

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
(5) MIN_MAPQ: "$MIN_MAPQ"
(6) SAMTOOLS_SIF: "$SAMTOOLS_SIF"
(7) THIS_ANALYSIS_DIR: "$THIS_ANALYSIS_DIR"
samtools version: "$SAMTOOLS_VERSION"
samtools view options:
-bh
-F 0x04
-q "$MIN_MAPQ"    # Is this too low?
"

EXPECTED_ARGS=7
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo "Usage: "$SCRIPT_TITLE" <SAMPLE_DIR> <ALIGNMENT_DIRNAME> <BAM_IN_FILENAME> <BAM_OUT_FILENAME> <MIN_MAPQ> <SAMTOOLS_SIF> <THIS_ANALYSIS_DIR>
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

BAM_OUT_FILE="$SAMPLE_DIR"/"$ALIGNMENT_DIRNAME"/"$BAM_OUT_FILENAME"

echo -e "${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]${NC}"

	srun singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view \
		-bh \
		-F 0x04 \
		-q "$MIN_MAPQ" \
		-o "$BAM_OUT_FILE" \
		"$BAM_IN_FILE"

	# check output status
	if [ $? -ne 0 ]
	then
			echo -e "${red}samtools view failed"
			exit 1
	fi

echo -e "${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]${NC}"

