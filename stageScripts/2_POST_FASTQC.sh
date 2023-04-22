#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=2_POST_FASTQC.sh
SCRIPT_VERSION=0.3
# DATE: 06-20-2018
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample (incl. both reads for PE) and is executed by <TTseq>_pipline.sh

# variables from command line via <TTseq>_pipline.sh:
SAMPLE_DIR=$1
SAMPLE_NAME=$2
FASTQ_FILE=$3
QC_DIR_NAME=$4
OUT_DIR_NAME=$5
THREADS=$6
# Other variables:
FASTQC_VERSION="$(fastqc -v | cut -d " " -f2)"

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
fastqc version: "$FASTQC_VERSION"
fastqc options:
-t 	"$THREADS"
--nogroup 	Disables grouping of bases for reads >50bp
--extract	zipped output file will be uncompressed
--outdir 	"$QC_DIR_NAME"/"$OUT_DIR_NAME"/"$SAMPLE_NAME"
"

EXPECTED_ARGS=6
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo -e "Usage: "$SCRIPT_TITLE" <SAMPLE_DIR> <SAMPLE_NAME> <FASTQ_FILE_NAME> <QC_DIR_NAME> <OUT_DIR_NAME> <THREADS>
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

cd "$QC_DIR_NAME"/"$OUT_DIR_NAME" # Simplified QC output by removing individual sample dirs

echo -e "
${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"

	srun fastqc -t "$THREADS" --nogroup --extract "$FASTQ_FILE" --outdir "$QC_DIR_NAME"/"$OUT_DIR_NAME" # Simplified QC output by removing individual sample dirs
	
	# check output status
	if [ $? -ne 0 ]
	then
			echo -e "${red}fastqc failed"
			exit 1
	fi

	# clean up - remove zip files
	rm $(basename ${FASTQ_FILE%.fastq.gz})_fastqc.zip
	
echo -e "${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]${NC}
"

