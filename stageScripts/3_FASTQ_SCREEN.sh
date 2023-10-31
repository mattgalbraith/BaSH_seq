#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=3_FASTQ_SCREEN.sh
SCRIPT_VERSION=0.4
# DATE: 06-21-2018
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample (incl. both reads for PE) and is executed by <XXseq>_pipline.sh
# Version 0.4 103023: updating to use FASTQC Singularity/Apptainer

# variables from command line via <TTseq>_pipline.sh:
SEQ_TYPE=${1}
SAMPLE_DIR=${2}
SAMPLE_NAME=${3}
FASTQR1_FILE=${4}
FASTQR2_FILE=${5}
QC_DIR_NAME=${6}
THREADS=${7}
FASTQSCREEN_CONF=${8}
FASTQSCREEN_SIF=${9}
THIS_ANALYSIS_DIR=${10}
RAW_DIR=${11}
REFS_DIR=${12}
# other variables:
OUT_DIR_NAME=fastq_screen
# FASTQ_SCREEN_CONF=$HOME/PipeLinesScripts/misc/fastq_screen.conf # now obtained as ARG 8
ALIGNER=bowtie2
# FASTQ_SCREEN_VERSION="$(fastq_screen --version | cut -d " " -f2)"
FASTQ_SCREEN_VERSION=`singularity run /data1/containers/fastqscreen0.15.2.sif bash -c 'fastq_screen --version  2> /dev/null'`

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for fastq_screen.sh:
(1) SEQ_TYPE: "$SEQ_TYPE"
(2) SAMPLE_DIR: "$SAMPLE_DIR"
(3) SAMPLE_NAME: "$SAMPLE_NAME"
(4) FASTQR1_FILE: "$FASTQR1_FILE"
(5) FASTQR2_FILE: "$FASTQR2_FILE"
(6) QC_DIR_NAME: "$QC_DIR_NAME"
(7) THREADS: "$THREADS"
(8) FASTQSCREEN_CONF: "$FASTQSCREEN_CONF"
(9) FASTQSCREEN_SIF: "$FASTQSCREEN_SIF"
(10) THIS_ANALYSIS_DIR: "$THIS_ANALYSIS_DIR"
(11) RAW_DIR: "$RAW_DIR"
(12) REFS_DIR: "$REFS_DIR"
fastq_screen version: "$FASTQ_SCREEN_VERSION"
fastq_screen options:
--conf "$FASTQSCREEN_CONF"
--aligner "$ALIGNER"
--threads "$THREADS"
--subset 1000000
--outdir
--paired
"

EXPECTED_ARGS=12
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo -e "Usage: "$SCRIPT_TITLE" <SEQ_TYPE> <SAMPLE_DIR> <SAMPLE_NAME> <FASTQR1_FILE> <FASTQR2_FILE> <QC_DIR_NAME> <THREADS> <FASTQSCREEN_CONF> <FASTQSCREEN_SIF> <THIS_ANALYSIS_DIR> <RAW_DIR> <REFS_DIR>
        ${red}ERROR - expecting "$EXPECTED_ARGS" ARGS but "$#" were provided:${NC}
        "$@"
        "
        exit 1
fi

# IF statements to check/create directories
	# check for read 1 FASTQ file
	if [ ! -f $FASTQR1_FILE ]
	then
		echo -e "${red}ERROR - file not found: "$FASTQR1_FILE"${NC}
		"
		exit 1
	fi

	# check if seq type is paired-end and read 2 FASTQ file
	if [ $SEQ_TYPE = "PE" ]
	then
		if [ ! -f $FASTQR2_FILE ]
		then
			echo -e  "${red}ERROR - file not found: "$FASTQR2_FILE"${NC}
			"
			echo
			exit 1
		fi
	fi

	# Check for fastq-screen output directory
	if [ ! -d "$QC_DIR_NAME"/fastq_screen/"$SAMPLE_NAME" ]
	then
		echo "making directory: "$QC_DIR_NAME"/fastq_screen/"$SAMPLE_NAME"
		"
		mkdir -p "$QC_DIR_NAME"/fastq_screen/"$SAMPLE_NAME"
	fi

cd $SAMPLE_DIR

echo -e "
${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"

	if [ "$SEQ_TYPE" = "SE" ]
		then
			echo "Running fastq_screen in single-end mode for "$SAMPLE_NAME"...
			"

			srun singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RAW_DIR":"$RAW_DIR" --bind "$REFS_DIR":"$REFS_DIR" --bind "$FASTQSCREEN_CONF":"$FASTQSCREEN_CONF" "$FASTQSCREEN_SIF" fastq_screen \
				--conf $FASTQSCREEN_CONF \
				--aligner $ALIGNER \
				--threads $THREADS \
				--subset 1000000 \
				--outdir "$QC_DIR_NAME"/fastq_screen/"$SAMPLE_NAME" \
				"$FASTQR1_FILE"
		
			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}fastq_screen single-end mode failed"
					exit 1
			fi
	elif [ "$SEQ_TYPE" = "PE" ]
		then
			echo "Running fastq_screen in paired-end mode for "$SAMPLE_NAME"...
			"

			srun singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RAW_DIR":"$RAW_DIR" --bind "$REFS_DIR":"$REFS_DIR" --bind "$FASTQSCREEN_CONF":"$FASTQSCREEN_CONF" "$FASTQSCREEN_SIF" fastq_screen \
				--conf $FASTQSCREEN_CONF \
				--aligner $ALIGNER \
				--threads $THREADS \
				--subset 1000000 \
				--outdir "$QC_DIR_NAME"/fastq_screen/"$SAMPLE_NAME" \
				--paired "$FASTQR1_FILE" "$FASTQR2_FILE"
		
			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}fastq_screen paired-end mode failed"
					exit 1
			fi

	else
		echo "${red}ERROR - SEQ_TYPE parameter not recognised:${NC} "$SEQ_TYPE"
		"
	fi

echo -e "
${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}
"

