#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=4_BWA-MEM.sh
SCRIPT_VERSION=0.1
# DATE: 11-02-2022
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample (incl. both reads for PE) and is executed by <CutRun>_pipeline.sh
# Reads are aligned using BWA MEM (running from Singularity container bwa-0.7.17-docker_amd64.sif) and directly converted to sorted BAM by samtools sort



# variables from command line via <CutRun>_pipeline.sh:
SEQ_TYPE=${1}
SAMPLE_DIR=${2}
SAMPLE_NAME=${3}
FASTQR1_FILE=${4}
FASTQR2_FILE=${5}
BWA_INDEX=${6}
BAM_OUT_FILENAME=${7}
OUT_DIR_NAME=${8}
THREADS=${9}
BWA_SIF=${10}
SAMTOOLS_SIF=${11}
THIS_ANALYSIS_DIR=${12}
RAW_DIR=${13}
REFS_DIR=${14}


# other variables
ALIGNMENT_SUMMARY_FILENAME="bwa-mem_alignment_summary.txt"
BWA_VERSION="$(singularity run "$BWA_SIF" bwa 2>&1 | grep Version)"
SAMTOOLS_VERSION=`singularity run "$SAMTOOLS_SIF" samtools --version | head -n1`

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for 4_BWA-MEM.sh:
(1) SEQ_TYPE: "$SEQ_TYPE"
(2) SAMPLE_DIR: "$SAMPLE_DIR"
(3) SAMPLE_NAME: "$SAMPLE_NAME"
(4) FASTQR1_FILE: "$FASTQR1_FILE"
(5) FASTQR2_FILE: "$FASTQR2_FILE"
(6) BWA_INDEX: "$BWA_INDEX"
(7) BAM_OUT_FILENAME: "$BAM_OUT_FILENAME"
(8) OUT_DIR_NAME: "$OUT_DIR_NAME"
(9) THREADS: "$THREADS"
(12) BWA_SIF: "$BWA_SIF"
(13) SAMTOOLS_SIF: "$SAMTOOLS_SIF"
(14) THIS_ANALYSIS_DIR: "$THIS_ANALYSIS_DIR"
(15) RAW_DIR: "$RAW_DIR"
(16) REFS_DIR: "$REFS_DIR"
bwa version: "$BWA_VERSION"
bwa options: mem
#
samtools version: "$SAMTOOLS_VERSION"
samtools sort options for BAM output:
-@ "$THREADS"
-m 32G
-T "$TMPDIR"
-o output file
samtools view options for counting number of aligned reads:
-F 256	do not output alignments with not primary flag 0x100
-F 4	do not output unmapped reads 0x4
"

# See: https://bio-bwa.sourceforge.net/bwa.shtml

EXPECTED_ARGS=14
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo -e "Usage: "$SCRIPT_TITLE" <SEQ_TYPE> <SAMPLE_DIR> <SAMPLE_NAME> <FASTQR1_FILE> <FASTQR2_FILE> <ALIGNER_INDEX> <BAM_OUT_FILENAME> <OUT_DIR_NAME> <THREADS> <BWA_SIF> <SAMTOOLS_SIF> <THIS_ANALYSIS_DIR> <RAW_DIR> <REFS_DIR>
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

	# Check for sample alignment output directory
	if [ ! -d $SAMPLE_DIR/$OUT_DIR_NAME ]
	then
		mkdir $SAMPLE_DIR/$OUT_DIR_NAME
	fi

	# Make temp directory for samtools
	SAMTOOLS_TMPDIR="$TMPDIR"/samtools/"$SAMPLE_NAME" # v0.3: updated for Proton2
	#
	if [ -d "$SAMTOOLS_TMPDIR" ]
	then
		rm -R "$SAMTOOLS_TMPDIR" 	# Better to overwrite in case error exit leaves previous version
	fi
	mkdir --parents "$SAMTOOLS_TMPDIR"



echo -e "${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"

	if [ "$SEQ_TYPE" = "SE" ]
		then
			echo "Running bwa-mem in single-end mode for "$SAMPLE_NAME"...
			"

			# () are used to group the bowtie2 and samtools commands to run in a subshell, allowing all output to be redirected to single file
			# set -o pipefail sets exit status of pipe/subshell to exit code of last command to exit non-zero (otherwise non-zero from bowtie2 will not be caught)
			# outer brackets with pipefail ensures that we do not just capture the exit code from tee
			(set -o pipefail
				(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RAW_DIR":"$RAW_DIR" --bind "$REFS_DIR":"$REFS_DIR" "$BWA_SIF" bwa mem \
				-t $THREADS \
				$BWA_INDEX \
				$FASTQR1_FILE | \
				singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools sort \
				-@ 4 -m 32G -T "$TMPDIR" -o "$SAMPLE_DIR/$OUT_DIR_NAME/$BAM_OUT_FILENAME" -) 2>&1 | tee "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$ALIGNMENT_SUMMARY_FILENAME")
				# The "-" at end of samtools view options tells it to use stdin as source file (piped directly from bowtie2)
				# OLD CMD: samtools view -bhS -o $SAMPLE_DIR/$OUT_DIR_NAME/$BAM_OUT_FILENAME -
				# samtools sort can convert directly to bam and removes the need for stage 7_SORT_BAM which uses slow Picard tool
			
			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}bwa mem or samtools failed in single-end mode"
					exit 1
			fi

	elif [ "$SEQ_TYPE" = "PE" ]
		then
			echo "Running hisat2 in paired-end mode for "$SAMPLE_NAME"...
			"
			# () are used to group the bowtie2 and samtools commands to run in a subshell, allowing all output to be redirected to single file
			# set -o pipefail sets exit status to exit code of pipe/subshell to last command to exit non-zero (otherwise non-zero from hisat2 will not be caught)
			# outer brackets with pipefail ensures that we do not just capture the exit code from tee
			(set -o pipefail
				(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RAW_DIR":"$RAW_DIR" --bind "$REFS_DIR":"$REFS_DIR" "$BWA_SIF" bwa mem \
				-t $THREADS \
				$BWA_INDEX \
				$FASTQR1_FILE \
				$FASTQR2_FILE | \
				singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools sort \
				-@ 4 -m 32G -T "$TMPDIR" -o "$SAMPLE_DIR/$OUT_DIR_NAME/$BAM_OUT_FILENAME" -) 2>&1 | tee "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$ALIGNMENT_SUMMARY_FILENAME")
			
			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}bwa mem or samtools failed in paired-end mode"
					exit 1
			fi

	else
			echo -e "${red}ERROR - SEQ_TYPE parameter not recognised:${NC} "$SEQ_TYPE"
			"
	fi

# samtools flagstat "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME" > "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME".flagstat.txt
# # count only primary alignments, this is the real mapped reads, flagstat is wrong?
# samtools view -F 256 -F 4 "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME" | wc -l > "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME".num_aligned.txt
# # these do not collect extra info / metrics for PE data - see PE_ALIGNMENT_METRICS.sh stage

echo -e "
${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}
"


