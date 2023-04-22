#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

# # Temporary export of environment variables to override R version until Bioconductor tools can be installed for R 3.6.1
# export RVERSION="3.3.3"
# export PATH=/gpfs/jespinosa/shared/Tools/R/$RVERSION/lib64/R/bin:$PATH
# export R_LIBS_SITE=$SHARED/Tools/lib/R_$RVERSION/
# export R_LIBS_USER=$SHARED/Tools/lib/R_$RVERSION/

SCRIPT_TITLE=11_HTSEQ_COUNT.sh
SCRIPT_VERSION=0.8
# DATE: 07-02-2021
# AUTHOR: Matthew Galbraith | Kohl Kinning
# SUMMARY: 
# This script is designed to be run once for each sample and is executed by <RNAseq>_pipeline.sh
# Version 0.3 removing srun command at line 122
# Version 0.4: adding settings for spike-in analysis
# Version 0.5: changing to sort by name for PE
# Version 0.6: lowering -m to 6G and changing over to new RPKM/TPM script
# Version 0.7: ?
# Version 0.8: modifications for Eureka: $TMPDIR -> /tmp; also had to add index (due to newer version of HTSeq-count?)

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
SEQ_TYPE=${15}
# other variables
BAMINFILE="$SAMPLE_DIR"/Alignments/"$SAMPLE_NAME".mapped.rgid.filtered.sorted.dups_mark.bam
HTSEQ_COUNT_VERSION="$(htseq-count -h | grep -e "framework")"
SAMTOOLS_VERSION="$(samtools --version 2>&1)"

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
(15) SEQ TYPE: "$SEQ_TYPE"
htseq-count version: "$HTSEQ_COUNT_VERSION"
htseq-count options:
--format=bam
--order="$SORT_ORDER"
--stranded=<yes/no/reverse>
--minaqual="$MINAQUAL"
--type="$FEATURETYPE"
--idattr="$IDATTR"
--mode="$MODE"
samtools version: "$SAMTOOLS_VERSION"
BAMINFILE: "$BAMINFILE"
"

EXPECTED_ARGS=15
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo -e "Usage: "$SCRIPT_TITLE" <ANALYSIS_DIR> <SAMPLE_DIR> <SAMPLE_NAME> <STRAND_TYPE> <SPIKE_IN> <SORT_ORDER> <GTF> <DM_GTF> <ERCC_GTF> <FEATURETYPE> <IDATTR> <MINAQUAL> <MODE> <COUNTS_DIRNAME> <SEQ_TYPE>
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
	COUNTS_DIR="$COUNTS_DIRNAME"/HTSeq
	if [ ! -d "$COUNTS_DIR" ]
	then
		echo "Creating directory: "$COUNTS_DIR""
		mkdir --parents "$COUNTS_DIR"
	fi
	OUTPUT_FILE="$COUNTS_DIR"/"$SAMPLE_NAME"_HTSeq_counts.txt

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

	# # Sort by name for PE to prevent buffer error in htseq-count
	# if [ "$SEQ_TYPE" = "PE" ] 
	# then
	# 	echo "Sorting bam file by name..."

	# 	SORT_ORDER="name" # used by HTseq count

	# 	samtools sort \
	# 		-m 6G \
	# 	 	-n \
	# 	 	-o /tmp/"$SAMPLE_NAME"_temp.bam \
	# 	 	-@ 16 \
	# 	 	"$BAMINFILE"

	# 	# check output status
	# 		if [ $? -ne 0 ]
	# 		then
	# 				echo -e "${red}samtools sort failed"
	# 				exit 1
	# 		else
	# 				echo "Done."
	# 		fi

# # Does not work with name sorted BAM
# 		samtools index \
# 		 	/tmp/"$SAMPLE_NAME"_temp.bam

	# 	# check output status
	# 		if [ $? -ne 0 ]
	# 		then
	# 				echo -e "${red}samtools index failed"
	# 				exit 1
	# 		else
	# 				echo "Done."
	# 		fi

	# 	BAMINFILE=/tmp/"$SAMPLE_NAME"_temp.bam

	# fi

	echo "Running htseq-count for main GTF..."

	htseq-count \
			--format=bam \
			--order="$SORT_ORDER" \
			--stranded="$STRANDED" \
			--minaqual="$MINAQUAL" \
			--type="$FEATURETYPE" \
			--idattr="$IDATTR" \
			--mode="$MODE" \
			"$BAMINFILE" \
			"$GTF" > ${OUTPUT_FILE%.txt}.tmp

			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}htseq-count failed"
					exit 1
			else
					echo "Done."
			fi

			# Extract special counters
			cat ${OUTPUT_FILE%.txt}.tmp | grep "^__" > ${OUTPUT_FILE%_counts.txt}_stats.txt
			# Remove empty gene names if they exist and special counters
			cat ${OUTPUT_FILE%.txt}.tmp | grep -v '^\s' | grep -v "^__" > "$OUTPUT_FILE"
			# Get total reads kept by HTSeq
			TOTAL_READS=$(awk '{s+=$2}END{print s}' "$OUTPUT_FILE")
			echo -e "Total_reads_kept:\t"$TOTAL_READS"" >> ${OUTPUT_FILE%_counts.txt}_stats.txt
			# mkdir "$ANALYSIS_DIR"/QC/HTSeq_count/
			# cp ${OUTPUT_FILE%_counts.txt}_stats.txt "$ANALYSIS_DIR"/QC/HTSeq_count/
			# Remove tmp file
			rm -rf ${OUTPUT_FILE%.txt}.tmp

		# TEMPORARILY INACTIVATED PENDING R SETUP ON EUREKA + UPDATING OF R SCRIPT
		# # generate TPMs
		# echo "Generating RPKMs and TPMs from counts file: "$OUTPUT_FILE""
		# srun Rscript $SHARED/Kohl/PipeLinesScripts/RNAseq/AnnotateCountsToRPKMAndTPM.R "$GTF" "$OUTPUT_FILE"  # NEED TO CHECK

		# check output status
		if [ $? -ne 0 ]
		then
				echo -e "${red}TPM file generation failed"
				exit 1
		fi

		# Count Drospohila reads if needed
		if [[ $SPIKE_IN = 'both' || $SPIKE_IN = 'dm' ]]
		then
			echo "Running htseq-count for Drosophila dm6 GTF..."

			htseq-count \
			--format=bam \
			--order="$SORT_ORDER" \
			--stranded="$STRANDED" \
			--minaqual="$MINAQUAL" \
			--type="$FEATURETYPE" \
			--idattr="$IDATTR" \
			--mode="$MODE" \
			"$BAMINFILE" \
			"$DM_GTF" > ${OUTPUT_FILE%.txt}_DM.tmp

			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}htseq-count failed"
					exit 1
			else
					echo "Done."
			fi

			# Extract special counters
			cat ${OUTPUT_FILE%.txt}_DM.tmp | grep "^__" > ${OUTPUT_FILE%_counts.txt}_DM_stats.txt
			# Remove empty gene names if they exist and special counters
			cat ${OUTPUT_FILE%.txt}_DM.tmp | grep -v '^\s' | grep -v "^__" > ${OUTPUT_FILE%.txt}_DM.txt
			# Get total reads kept by HTSeq
			TOTAL_READS=$(awk '{s+=$2}END{print s}' ${OUTPUT_FILE%.txt}_DM.txt)
			echo -e "Total_reads_kept:\t"$TOTAL_READS"" >> ${OUTPUT_FILE%_counts.txt}_DM_stats.txt
			# Remove tmp file
			rm -rf ${OUTPUT_FILE%.txt}_DM.tmp

		fi

		# Count ERCC reads if needed
		if [[ $SPIKE_IN = 'both' || $SPIKE_IN = 'ercc' ]]
		then
			echo "Running htseq-count for ERCC GTF..."

			htseq-count \
			--format=bam \
			--order="$SORT_ORDER" \
			--stranded="$STRANDED" \
			--minaqual="$MINAQUAL" \
			--type="$FEATURETYPE" \
			--idattr="$IDATTR" \
			--mode="$MODE" \
			"$BAMINFILE" \
			"$ERCC_GTF" > ${OUTPUT_FILE%.txt}_ERCC.tmp

			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}htseq-count failed"
					exit 1
			else
					echo "Done."
			fi

			# Extract special counters
			cat ${OUTPUT_FILE%.txt}_ERCC.tmp | grep "^__" > ${OUTPUT_FILE%_counts.txt}_ERCC_stats.txt
			# Remove empty gene names if they exist and special counters
			cat ${OUTPUT_FILE%.txt}_ERCC.tmp | grep -v '^\s' | grep -v "^__" > ${OUTPUT_FILE%.txt}_ERCC.txt
			# Get total reads kept by HTSeq
			TOTAL_READS=$(awk '{s+=$2}END{print s}' ${OUTPUT_FILE%.txt}_ERCC.txt)
			echo -e "Total_reads_kept:\t"$TOTAL_READS"" >> ${OUTPUT_FILE%_counts.txt}_ERCC_stats.txt
			# Remove tmp file
			rm -rf ${OUTPUT_FILE%.txt}_ERCC.tmp
			rm -rf "$BAMINFILE"

		fi


echo -e "
${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}
"


