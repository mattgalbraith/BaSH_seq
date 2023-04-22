#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=14_SALMON.sh
SCRIPT_VERSION=0.1
# DATE: 12-05-2022
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample (incl. both reads for PE) and is executed by <RNAseq>_pipeline.sh
# Reads are quantified at transcript level using Salmon (running from Singularity container).



# variables from command line via <CutRun>_pipeline.sh:
SEQ_TYPE=${1}
SAMPLE_DIR=${2}
SAMPLE_NAME=${3}
FASTQR1_FILE=${4}
FASTQR2_FILE=${5}
SALMON_INDEX=${6}
OUT_DIR_NAME=${7}
THREADS=${8}

# other variables
# ALIGNMENT_SUMMARY_FILENAME="bwa-mem_alignment_summary.txt"
SALMON_SIF=/home/matthew.galbraith/Tools/Singularity/salmon.sif
SALMON_VERSION="$(singularity run $SALMON_SIF salmon --version)"

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for 14_SALMON.sh:
(1) SEQ_TYPE: "$SEQ_TYPE"
(2) SAMPLE_DIR: "$SAMPLE_DIR"
(3) SAMPLE_NAME: "$SAMPLE_NAME"
(4) FASTQR1_FILE: "$FASTQR1_FILE"
(5) FASTQR2_FILE: "$FASTQR2_FILE"
(6) SALMON_INDEX: "$SALMON_INDEX"
(7) OUT_DIR_NAME: "$OUT_DIR_NAME"
(8) THREADS: "$THREADS"
Salmon version: "$SALMON_VERSION"
Salmon options:
see Salmon log files
"

# See: https://salmon.readthedocs.io/en/latest/

EXPECTED_ARGS=8
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo -e "Usage: "$SCRIPT_TITLE" <SEQ_TYPE> <SAMPLE_DIR> <SAMPLE_NAME> <FASTQR1_FILE> <FASTQR2_FILE> <SALMON_INDEX> <OUT_DIR_NAME> <THREADS>
        ${red}ERROR - expecting "$EXPECTED_ARGS" ARGS but "$#" were provided:${NC}
        "$@"
        "
        exit 1
fi

# IF statements to check/create directories
	# Check for Sample directory
	if [ ! -d "$SAMPLE_DIR" ]
	then
	        echo -e "${red}ERROR - sample directory not found:${NC} "$SAMPLE_DIR"
	        "
	        exit 1
	fi

	# check for read 1 FASTQ file
	if [ ! -f $FASTQR1_FILE ]
	then
		echo -e "${red}ERROR - File not found: "$FASTQR1_FILE"
		"
		exit 1
	fi

	# check if seq type is paired-end and read 2 FASTQ file
	if [ $SEQ_TYPE = "PE" ]
	then
		if [ ! -f $FASTQR2_FILE ]
		then
			echo -e "${red}ERROR - File not found: "$FASTQR2_FILE"
			"
			exit 1
		fi
	fi

	# Check for sample Salmon output directory
	if [ ! -d $SAMPLE_DIR/$OUT_DIR_NAME ]
	then
		mkdir $SAMPLE_DIR/$OUT_DIR_NAME
	fi



echo -e "${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"

	if [ "$SEQ_TYPE" = "SE" ]
		then
			echo "Running Salmon in single-end mode for "$SAMPLE_NAME"...
			"

			echo "Single-end mode not yet implemented for Salmon"

			exit 1
			
			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}Salmon failed in single-end mode"
					exit 1
			fi

	elif [ "$SEQ_TYPE" = "PE" ]
		then
			echo "Running Salmon in paired-end mode for "$SAMPLE_NAME"...
			"
			
			singularity run $SALMON_SIF salmon quant \
        			-p $THREADS \
        			-i $SALMON_INDEX \
        			--libType A \
        			--numGibbsSamples 30 \
        			--seqBias \
        			--gcBias \
        			-o "$OUT_DIR_NAME" \
        			-1 $FASTQR1_FILE \
        			-2 $FASTQR2_FILE
			
			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}Salmon failed in paired-end mode"
					exit 1
			fi

	else
			echo -e "${red}ERROR - SEQ_TYPE parameter not recognised:${NC} "$SEQ_TYPE"
			"
	fi


echo -e "
${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}
"


