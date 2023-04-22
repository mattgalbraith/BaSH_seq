#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=11_FEATURECOUNTS_RNA_v0.1.sh
SCRIPT_VERSION=0.1
# DATE: 05-10-2019
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample and is executed by <RNAseq>_pipeline.sh
## CHANGE LOG ##

# variables from command line via <RNAseq>_pipeline.sh:
ANALYSIS_DIR=${1}
SAMPLE_DIR=${2}
SAMPLE_NAME=${3}
STRAND_TYPE=${4}
SPIKE_IN=${5}
SORT_ORDER=${6}
GTF=${7}
DM_GTF=${8}
ERCC_GTF=$9
FEATURETYPE=${10}
IDATTR=${11}
MINAQUAL=${12}
MODE=${13}
COUNTS_DIRNAME=${14}
# other variables
BAMINFILE="$SAMPLE_DIR"/Alignments/"$SAMPLE_NAME".mapped.rgid.filtered.sorted.dups_mark.bam
FEATURECOUNTS_VERSION="$(featureCounts -v 2>&1 | cut -d " " -f2)"

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for "$SCRIPT_TITLE":
(1) ANALYSIS_DIR: "$ANALYSIS_DIR"
(2) SAMPLE_DIR: "$SAMPLE_DIR"
(3) SAMPLE_NAME: "$SAMPLE_NAME"
(4) STRAND_TYPE: "$STRAND_TYPE"
(5) SPIKE_IN: "$SPIKE_IN"
(6) SORT_ORDER: "$SORT_ORDER"
(7) GTF: "$GTF"
(8) Drosophila dm6 GTF: "$DM_GTF"
(9) ERCC GTF: "$ERCC_GTF"
(10) FEATURETYPE: "$FEATURETYPE"
(11) IDATTR: "$IDATTR"
(12) MINAQUAL: "$MIN_QUAL"
(13) MODE: "$MODE"
(14) COUNTS_DIRNAME: "$COUNTS_DIRNAME"
featureCounts version: "$FEATURECOUNTS_VERSION"
featureCounts options:
-a "$GTF"
-F GTF
-t "$FEATURETYPE"
-g "$IDATTR"
-Q "$MINAQUAL"
-p
-s 2
-C
--primary
--tmpDir $TMPDIR
-T 12
--verbose
"

EXPECTED_ARGS=14
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo -e "Usage: "$SCRIPT_TITLE" <ANALYSIS_DIR> <SAMPLE_DIR> <SAMPLE_NAME> <STRAND_TYPE> <SPIKE_IN> <SORT_ORDER> <GTF> <DM_GTF> <ERCC_GTF> <FEATURETYPE> <IDATTR> <MINAQUAL> <MODE> <COUNTS_DIRNAME>
        ${red}ERROR - expecting "$EXPECTED_ARGS" ARGS but "$#" were provided:${NC}
        "$@"
        "
        exit 1
fi

# IF statements to check/create directories
	# Check for Sample Alignment directory
	if [ ! -d "$SAMPLE_DIR"/Alignments ]
	then
		echo "${red}ERROR - folder does not exist: "$SAMPLE_DIR"/Alignments
		"
		exit 1
	
	elif [ -d "$SAMPLE_DIR"/Alignments ]
	then
		echo "Alignments folder found: "$SAMPLE_DIR"/Alignments
		"
	fi

	# Check for sample counts output directory
	COUNTS_DIR="$COUNTS_DIRNAME"/featureCounts
	if [ ! -d "$COUNTS_DIR" ]
	then
		echo "Creating directory: "$COUNTS_DIR""
		mkdir --parents "$COUNTS_DIR"
	fi
	OUTPUT_FILE="$COUNTS_DIR"/"$SAMPLE_NAME"_featureCounts.txt

# Set stranded option
	# This ensures that reads are assigned correctly to fwd/rev strand genes
	# Check RSeQC output to confirm correct setting
	if [ "$STRAND_TYPE" = "unstranded" ]
	then
		echo "Setting --stranded to "no""
		STRANDED="no"
	elif [ "$STRAND_TYPE" = "strand-specific-fwd" ]
		then
			echo "Setting --stranded to "yes"" # eg NuGen Universal RNA-seq libraries CHECK
			STRANDED="yes"
	elif [ "$STRAND_TYPE" = "strand-specific-rev" ]
		then
			echo "Setting --stranded to "reverse"" # eg Illumina TruSeq RNA libraries CHECK
			STRANDED="reverse"
	else
		echo "${red}ERROR - Strand type not recognized:${NC} "$STRAND_TYPE"
		"
		exit 1
	fi



echo -e "
${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"

	echo "Running featureCounts for main GTF..."

	# note -p only for PE
	# note -s only for stranded
	# note could also set -B to require both ends mapped

	srun featureCounts \
			-a "$GTF" \
			-F GTF \
			-t "$FEATURETYPE" \
			-g "$IDATTR" \
			-Q "$MINAQUAL" \
			-p \
			-s 2 \
			-C \
			--primary \
			--tmpDir $TMPDIR \
			-T 12 \
			--verbose \
			-o "$OUTPUT_FILE" \
			"$BAMINFILE"

		# check output status
		if [ $? -ne 0 ]
		then
				echo -e "${red}featureCounts failed"
				exit 1
		else
				echo "Done."
		fi


		# TURNED OFF FOR NOW - need to skip first line for featureCounts output or do later in R
		# # generate RPKMs
		# echo "Generating RPKMs from counts file: "$OUTPUT_FILE""
		# srun Rscript $SHARED/Matt/PipeLinesScripts/RNAseq/AnnotateCountsToRPKMs_v2.R "$GTF" "$OUTPUT_FILE"

		# check output status
		if [ $? -ne 0 ]
		then
				echo -e "${red}RPKM file generation failed"
				exit 1
		fi

		# Count Drospohila reads if needed
		if [[ $SPIKE_IN = 'both' || $SPIKE_IN = 'dm' ]]
		then
			echo "Running featureCounts for Drosophila dm6 GTF..."

			srun featureCounts \
				-a "$DM_GTF" \
				-F GTF \
				-t "$FEATURETYPE" \
				-g "$IDATTR" \
				-Q "$MINAQUAL" \
				-p \
				-s 2 \
				-C \
				--primary \
				--tmpDir $TMPDIR \
				-T 16 \
				--verbose \
				-o ${OUTPUT_FILE%.txt}_DM.txt \
				"$BAMINFILE"

			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}featureCounts failed"
					exit 1
			else
					echo "Done."
			fi

		fi

		# Count ERCC reads if needed
		if [[ $SPIKE_IN = 'both' || $SPIKE_IN = 'ercc' ]]
		then
			echo "Running featureCounts for ERCC GTF..."

			srun featureCounts \
				-a "$ERCC_GTF" \
				-F GTF \
				-t "$FEATURETYPE" \
				-g "$IDATTR" \
				-Q "$MINAQUAL" \
				-p \
				-s 2 \
				-C \
				--primary \
				--tmpDir $TMPDIR \
				-T 12 \
				--verbose \
				-o ${OUTPUT_FILE%.txt}_ERCC.txt \
				"$BAMINFILE"

			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}featureCounts failed"
					exit 1
			else
					echo "Done."
			fi

		fi


echo -e "
${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}
"

