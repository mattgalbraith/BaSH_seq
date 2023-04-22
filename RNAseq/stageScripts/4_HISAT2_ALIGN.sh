#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=4_HISAT2-DNA.sh
SCRIPT_VERSION=0.2
# DATE: 01-29-2020
# AUTHOR: Matthew Galbraith, Kohl Kinning
# SUMMARY: 
# This script is designed to be run once for each sample (incl. both reads for PE) and is executed by <TTseq>_pipline.sh
# Reads are aligned using HISAT2 and directly converted to sorted BAM by samtools sort
# Version 0.2: Removed --no-spliced-alignment and --maxins parameters for general RNASeq use. Added splice-sites param, a file generated when building the index



# variables from command line via <TTseq>_pipeline.sh:
SEQ_TYPE=${1}
STRAND_TYPE=${2}
SAMPLE_DIR=${3}
SAMPLE_NAME=${4}
FASTQR1_FILE=${5}
FASTQR2_FILE=${6}
HISAT2_INDEX=${7}
HISAT2_SS=${8}
BAM_OUT_FILENAME=${9}
OUT_DIR_NAME=${10}
THREADS=${11}

# STRAND_TYPE=${9}
# MATE_INNER_DIST=${10}

# other variables
ALIGNMENT_SUMMARY_FILENAME="hisat_alignment_summary.txt"
HISAT2_VERSION="$(hisat2 --version | head -n1 | awk '{print $3}')"
SAMTOOLS_VERSION="$(samtools --version 2>&1)"

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for 4_HISAT2-DNA.sh:
(1) SEQ_TYPE: "$SEQ_TYPE"
(2) STRAND_TYPE: "$STRAND_TYPE"
(3) SAMPLE_DIR: "$SAMPLE_DIR"
(4) SAMPLE_NAME: "$SAMPLE_NAME"
(5) FASTQR1_FILE: "$FASTQR1_FILE"
(6) FASTQR2_FILE: "$FASTQR2_FILE"
(7) HISAT2_INDEX: "$HISAT2_INDEX"
(8) HISAT2_SS: "$HISAT2_SS"
(9) BAM_OUT_FILENAME: "$BAM_OUT_FILENAME"
(10) OUT_DIR_NAME: "$OUT_DIR_NAME"
(11) THREADS: "$THREADS"
hisat2 version: "$HISAT2_VERSION"
hisat2 options:
# --phred33
# --sensitive = default
# --end-to-end = default
# -x 		basename of the index for the reference genome
# -U 		unpaired (single-end) input fastq file
# -1		read 1 (paired-end) input fastq file
# -2		read 2 (paired-end) input fastq file
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

# See: https://ccb.jhu.edu/software/hisat2/manual.shtml

EXPECTED_ARGS=11
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo -e "Usage: "$SCRIPT_TITLE" <SEQ_TYPE> <STRAND_TYPE> <SAMPLE_DIR> <SAMPLE_NAME> <FASTQR1_FILE> <FASTQR2_FILE> <ALIGNER_INDEX> <BAM_OUT_FILENAME> <OUT_DIR_NAME> <THREADS>
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

	# Set mate orientation option: (this should not matter for TT-seq until later steps)
	# The upstream/downstream mate orientations for a valid paired-end alignment against the forward reference strand. E.g., if --fr is specified and there is a candidate paired-end alignment where mate 1 appears upstream of the reverse complement of mate 2 and the fragment length constraints (-I and -X) are met, that alignment is valid. Also, if mate 2 appears upstream of the reverse complement of mate 1 and all other constraints are met, that too is valid. --rf likewise requires that an upstream mate1 be reverse-complemented and a downstream mate2 be forward-oriented. --ff requires both an upstream mate 1 and a downstream mate 2 to be forward-oriented. Default: --fr (appropriate for Illumina's Paired-end Sequencing Assay).
	if [ "$STRAND_TYPE" = "unstranded" ]
		then
			LIB_TYPE=""
	elif [ "$STRAND_TYPE" = "strand-specific-fwd" ]
		then
			echo "Setting mate orientation to "--rf"" # eg NuGen Universal PE libraries CHECK
			LIB_TYPE="--rf"
	elif [ "$STRAND_TYPE" = "strand-specific-rev" ]
		then
			echo "Setting --mate orientation to "--fr"" # eg Illumina PE libraries
			LIB_TYPE="--fr"
	else
		echo "${red}ERROR - Strand type not recognized:${NC} "$STRAND_TYPE"
		"
		exit 1
	fi


echo -e "${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"

	if [ "$SEQ_TYPE" = "SE" ]
		then
			echo "Running hisat2 in single-end mode for "$SAMPLE_NAME"...
			"

			# () are used to group the bowtie2 and samtools commands to run in a subshell, allowing all output to be redirected to single file
			# set -o pipefail sets exit status of pipe/subshell to exit code of last command to exit non-zero (otherwise non-zero from bowtie2 will not be caught)
			# outer brackets with pipefail ensures that we do not just capture the exit code from tee
			(set -o pipefail
				(hisat2 \
				-p $THREADS \
				-x "$HISAT2_INDEX" \
 				--known-splicesite-infile "$HISAT2_SS" \
				--phred33 \
				-U $FASTQR1_FILE | \
				samtools sort -@ 12 -m 32G -T "$TMPDIR" -o "$SAMPLE_DIR/$OUT_DIR_NAME/$BAM_OUT_FILENAME" -) 2>&1 | tee "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$ALIGNMENT_SUMMARY_FILENAME")
				# The "-" at end of samtools view options tells it to use stdin as source file (piped directly from bowtie2)
				# OLD CMD: samtools view -bhS -o $SAMPLE_DIR/$OUT_DIR_NAME/$BAM_OUT_FILENAME -
				# samtools sort can convert directly to bam and removes the need for stage 7_SORT_BAM which uses slow Picard tool
			
			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}hisat2 or samtools failed in single-end mode"
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
				(hisat2 \
				-p $THREADS \
				-x "$HISAT2_INDEX" \
				--known-splicesite-infile "$HISAT2_SS" \
				--phred33 \
				"$LIB_TYPE" \
				-1 $FASTQR1_FILE \
				-2 $FASTQR2_FILE | \
				samtools sort -@ 12 -m 32G -T "$TMPDIR" -o "$SAMPLE_DIR/$OUT_DIR_NAME/$BAM_OUT_FILENAME" -) 2>&1 | tee "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$ALIGNMENT_SUMMARY_FILENAME")
			
			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}hisat2 or samtools failed in paired-end mode"
					exit 1
			fi

	else
			echo -e "${red}ERROR - SEQ_TYPE parameter not recognised:${NC} "$SEQ_TYPE"
			"
	fi

samtools flagstat "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME" > "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME".flagstat.txt
# count only primary alignments, this is the real mapped reads, flagstat is wrong?
samtools view -F 256 -F 4 "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME" | wc -l > "$SAMPLE_DIR"/"$OUT_DIR_NAME"/"$BAM_OUT_FILENAME".num_aligned.txt
# these do not collect extra info / metrics for PE data - see PE_ALIGNMENT_METRICS.sh stage

echo -e "
${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}
"


