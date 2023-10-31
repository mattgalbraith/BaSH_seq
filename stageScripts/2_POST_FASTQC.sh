#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=2_POST_FASTQC.sh
SCRIPT_VERSION=0.4
# DATE: 06-20-2018
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample (incl. both reads for PE) and is executed by <TTseq>_pipline.sh
# Version 0.4 102623: updating to use FASTQC Singularity/Apptainer

# variables from command line via <TTseq>_pipline.sh:
SAMPLE_DIR=$1
SAMPLE_NAME=$2
FASTQ_FILE=$3
QC_DIR_NAME=$4
OUT_DIR_NAME=$5
THREADS=$6
FASTQC_SIF=$7
THIS_ANALYSIS_DIR=$8
RAW_DIR=$9
# Other variables:
# FASTQC_VERSION="$(fastqc -v | cut -d " " -f2)"
FASTQC_VERSION=`singularity run "$FASTQC_SIF" bash -c 'fastqc --version 2>NULL' | cut -d " " -f2,2`


blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for FASTQC.sh:
(1) SAMPLE_DIR: "$SAMPLE_DIR"
(2) SAMPLE_NAME: "$SAMPLE_NAME"
(3) FASTQ_FILE: "$FASTQ_FILE"
(4) QC_DIR_NAME: "$QC_DIR_NAME"
(5) OUT_DIR_NAME: "$OUT_DIR_NAME"
(6) THREADS: "$THREADS" (-t option)
(7) FASTQC container: "$FASTQC_SIF"
(8) THIS_ANALYSIS_DIR: "$THIS_ANALYSIS_DIR"
(9) RAW_DIR: "$RAW_DIR"
fastqc version: "$FASTQC_VERSION"
fastqc options:
-t 	"$THREADS"
--nogroup 	Disables grouping of bases for reads >50bp
--extract	zipped output file will be uncompressed
--outdir 	"$QC_DIR_NAME"/"$OUT_DIR_NAME"/"$SAMPLE_NAME"
"

EXPECTED_ARGS=9
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo -e "Usage: "$SCRIPT_TITLE" <SAMPLE_DIR> <SAMPLE_NAME> <FASTQ_FILE_NAME> <QC_DIR_NAME> <OUT_DIR_NAME> <THREADS> <FASTQC_SIF> <THIS_ANALYSIS_DIR> <RAW_DIR>
        ${red}ERROR - expecting "$EXPECTED_ARGS" ARGS but "$#" were provided:${NC}
        "$@"
        "
        exit 1
fi

# IF statements to check/create directories
	# Check for Sample directory
	if [ ! -d "$SAMPLE_DIR" ]
	then
	        echo -e "${red}ERROR - sample directory not found:${NC} "$SAMPLE_DIR""
	        exit 1
	fi

	# Check for FASTQC output directory
	if [ ! -d "$QC_DIR_NAME"/"$OUT_DIR_NAME" ] # Simplified QC output by removing individual sample dirs
	then
		echo "making directory: "$QC_DIR_NAME"/"$OUT_DIR_NAME"" # Simplified QC output by removing individual sample dirs
		mkdir --parents "$QC_DIR_NAME"/"$OUT_DIR_NAME" # Simplified QC output by removing individual sample dirs
	fi

# cd "$QC_DIR_NAME"/"$OUT_DIR_NAME" # Simplified QC output by removing individual sample dirs # NOT SURE WHY CHANGING TO THIS DIR
 # + prevents Apptainer auto mounting Analysis dir as PWD
cd "$THIS_ANALYSIS_DIR" # THIS DOES NOT APPEAR TO BE SUFFICENT TO AUTO MOUNT ANALYSIS_DIR

echo -e "
${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"

	srun singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RAW_DIR":"$RAW_DIR" "$FASTQC_SIF" fastqc -t "$THREADS" --nogroup --extract "$FASTQ_FILE" --outdir "$QC_DIR_NAME"/"$OUT_DIR_NAME"		# Simplified QC output by removing individual sample dirs
	# Need to mount (bind) RAW_DIR so that apptainer can read symbolic links to FASTQ files
	
	# check output status
	if [ $? -ne 0 ]
	then
			echo -e "${red}fastqc failed"
			exit 1
	fi

	# clean up - remove zip files
	rm "$QC_DIR_NAME"/"$OUT_DIR_NAME"/$(basename ${FASTQ_FILE%.fastq.gz})_fastqc.zip
	
echo -e "${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]${NC}
"

