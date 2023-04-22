#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=4_TOPHAT2_ALIGN.sh
SCRIPT_VERSION=0.3
# DATE: 09-17-2019
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample (incl. both reads for PE) and is executed by <RNAseq>_pipeline.sh
# Reads are aligned to reference genome and exon annotation using Tophat2
#
# Change log
# v0.3: changing ln -s to ln -sf to force link creation if already existing

# variables from command line via <RNAseq>_pipeline.sh:
SEQ_TYPE=${1}
SAMPLE_DIR=${2}
SAMPLE_NAME=${3}
FASTQR1_FILE=${4}
FASTQR2_FILE=${5}
BOWTIE_INDEX=${6}
GTF=${7}
TRANSCRIPTOME_INDEX=${8}
STRAND_TYPE=${9}
MATE_INNER_DIST=${10}
MATE_STD_DEV=${11}
BAM_OUT_FILENAME=${12}
OUT_DIR_NAME=${13}
THREADS=${14}
# other variables
# ALIGNMENT_SUMMARY_FILENAME="bowtie_alignment_summary.txt"
TOPHAT2_VERSION="echo $(tophat2 --version | cut -d " " -f2)"
SAMTOOLS_VERSION="$(samtools --version 2>&1)"

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for "$SCRIPT_TITLE":
(1) SEQ_TYPE: "$SEQ_TYPE"
(2) SAMPLE_DIR: "$SAMPLE_DIR"
(3) SAMPLE_NAME: "$SAMPLE_NAME"
(4) FASTQR1_FILE: "$FASTQR1_FILE"
(5) FASTQR2_FILE: "$FASTQR2_FILE"
(6) BOWTIE_INDEX: "$BOWTIE_INDEX"
(7) GTF: "$GTF"
(8) TRANSCRIPTOME_INDEX: "$TRANSCRIPTOME_INDEX"
(9) STRAND_TYPE: "$STRAND_TYPE"
(10) MATE_INNER_DIST: "$MATE_INNER_DIST"
(11) MATE_STD_DEV: "$MATE_STD_DEV"
(12) BAM_OUT_FILENAME: "$BAM_OUT_FILENAME"
(13) OUT_DIR_NAME: "$OUT_DIR_NAME"
(14) THREADS: "$THREADS"
tophat2 version: "$TOPHAT2_VERSION"
tophat2 options:
-p "$THREADS"
--b2-sensitive
--keep-fasta-order
--no-coverage-search
--max-multihits 10
--library-type "$TOPHAT_STRAND_OPTS"
--GTF "$GTF"
--transcriptome-index="$TRANSCRIPTOME_INDEX"
samtools version: "$SAMTOOLS_VERSION"
"

EXPECTED_ARGS=14
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo -e "Usage: "$SCRIPT_TITLE" <SEQ_TYPE> <SAMPLE_DIR> <SAMPLE_NAME> <FASTQR1_FILE> <FASTQR2_FILE> <BOWTIE_INDEX> <GTF> <TRANSCRIPTOME_INDEX> <STRAND_TYPE> <MATE_INNER_DIST> <MATE_STD_DEV> <BAM_OUT_FILENAME> <OUT_DIR_NAME> <THREADS>
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

	# Check for sample tophat2 alignment output directory
	if [ ! -d "$SAMPLE_DIR"/"$OUT_DIR_NAME" ]
	then
		echo "Creating directory: "$SAMPLE_DIR"/"$OUT_DIR_NAME""
		mkdir "$SAMPLE_DIR"/"$OUT_DIR_NAME"
	fi

	# Set library-type option
	# This ensures that the correct XS attribute tags are set for each read in SAM/BAM file (ie will not affect actual alignment but will affect downstream gene-based analysis)
	# Check RSeQC output to confirm correct setting
	if [ "$STRAND_TYPE" = "unstranded" ]
		then
			echo "Setting --library-type to "fr-unstranded""
			LIB_TYPE="fr-unstranded"
	elif [ "$STRAND_TYPE" = "strand-specific-fwd" ]
		then
			echo "Setting --library-type to "fr-secondstrand"" # eg NuGen Universal RNA-seq libraries CHECK
			LIB_TYPE="fr-secondstrand"
	elif [ "$STRAND_TYPE" = "strand-specific-rev" ]
		then
			echo "Setting --library-type to "fr-firststrand"" # eg Illumina TruSeq RNA libraries CHECK
			LIB_TYPE="fr-firststrand"
	else
		echo "${red}ERROR - Strand type not recognized:${NC} "$STRAND_TYPE"
		"
		exit 1
	fi


cd "$SAMPLE_DIR"/"$OUT_DIR_NAME"

echo -e "${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"

	if [ "$SEQ_TYPE" = "SE" ]
		then
			echo "Running tophat2 in single-end mode for "$SAMPLE_NAME"...
			"
	
			srun tophat2 -p "$THREADS" \
				--b2-sensitive \
				--keep-fasta-order \
				--no-coverage-search \
				--max-multihits 10 \
				--library-type "$LIB_TYPE" \
				--GTF "$GTF" \
				--transcriptome-index="$TRANSCRIPTOME_INDEX" \
				"$BOWTIE_INDEX" \
				"$FASTQR1_FILE"

		# check output status
		if [ $? -ne 0 ]
		then
				echo -e "${red}tophat2 failed in single-end mode${NC}"
				exit 1
		fi
	
		# Generate symbolic links
		ln -sf "$SAMPLE_DIR"/"$OUT_DIR_NAME"/tophat_out/accepted_hits.bam "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME"
		if [ $? -ne 0 ]
		then
			echo "${red}ERROR - Problem with BAM output file:${NC} "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME""
			exit 1
		fi
		ln -sf "$SAMPLE_DIR"/"$OUT_DIR_NAME"/tophat_out/align_summary.txt "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME".tophat2_alignment_summary.txt

	elif [ "$SEQ_TYPE" = "PE" ]
		then
			echo "Running tophat22 in paired-end mode for "$SAMPLE_NAME"...
			"
			# --mate-inner-dist and --mate-std-dev 20 should be set based on library information if possible but may have to be set after an initial alignment run
			srun tophat2 -p "$THREADS" \
				--b2-sensitive \
				--keep-fasta-order \
				--no-coverage-search \
				--max-multihits 10 \
				--library-type "$LIB_TYPE" \
				--GTF "$GTF" \
				--transcriptome-index="$TRANSCRIPTOME_INDEX" \
				--mate-inner-dist "$MATE_INNER_DIST" \
				--mate-std-dev "$MATE_STD_DEV" \
				"$BOWTIE_INDEX" \
				"$FASTQR1_FILE" "$FASTQR2_FILE"

			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}tophat2 failed in paired-end mode${NC}"
					exit 1
			fi

			# Generate symbolic links
			ln -sf "$SAMPLE_DIR"/"$OUT_DIR_NAME"/tophat_out/accepted_hits.bam "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME"
			if [ $? -ne 0 ]
			then
				echo "${red}ERROR - Problem with BAM output file:${NC} "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME""
				exit 1
			fi
			ln -sf "$SAMPLE_DIR"/"$OUT_DIR_NAME"/tophat_out/align_summary.txt "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME".tophat2_alignment_summary.txt

	else
			echo -e "${red}ERROR - SEQ_TYPE parameter not recognised:${NC} "$SEQ_TYPE"
			"
			exit 1

	fi

srun samtools flagstat "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME" > "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME".flagstat.txt
# count only primary alignments, this is the real mapped reads, flagstat is wrong?
srun samtools view -F 256 -F 4 "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME" | wc -l > "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME".num_aligned.txt
# these do not collect extra info / metrics for PE data - see PE_ALIGNMENT_METRICS.sh stage

echo -e "
${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}
"

