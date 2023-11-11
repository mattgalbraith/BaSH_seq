#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=10_RSEQC.sh
SCRIPT_VERSION=0.6
# DATE: 06-21-2018
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample and is executed by <RNAseq>_pipline.sh
#
# v0.4 070921
# changes to use reseqc v4.0.0 Singularity container made by John Finigan
# v0.5 090922
# correction for tin section to use Housekeeping BED file instead of RefSeq BED (line 310)
# Version 0.4 111023: 
# updating to use RSEQC Singularity/Apptainer on Proton2


# variables from command line via <RNAseq>_pipline.sh:
SEQ_TYPE=${1}
THIS_ANALYSIS_DIR=${2}
QC_DIR=${3}
SAMPLE_DIR=${4}
ALIGNMENT_DIRNAME=${5}
BAM_IN_FILENAME=${6}
SAMPLE_NAME=${7}
REFSEQ_BED=${8}
HOUSEKEEPING_BED=${9}
RSEQC_SIF=${10}
RSEQC_REFS=${11}
# other variables:
# RSEQC_VERSION="$(ls /gpfs/jespinosa/shared/Tools/omics_tools/rseqc/)"
RSEQC_VERSION="$(singularity run "$RSEQC_SIF" bam_stat.py --version | sed 's/bam_stat.py\s//')"
MIN_MAP_QUAL=10

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for "$SCRIPT_TITLE":
(1) SEQ_TYPE: "$SEQ_TYPE"
(2) THIS_ANALYSIS_DIR: "$THIS_ANALYSIS_DIR"
(3) QC_DIR: "$QC_DIR"
(4) SAMPLE_DIR: "$SAMPLE_DIR"
(5) ALIGNMENT_DIRNAME: "$ALIGNMENT_DIRNAME"
(6) BAM_IN_FILENAME: "$BAM_IN_FILENAME"
(7) SAMPLE_NAME: "$SAMPLE_NAME"
(8) REFSEQ_BED: "$REFSEQ_BED"
(9) HOUSEKEEPING_BED: "$HOUSEKEEPING_BED"
(10) RSEQC_SIF: "$RSEQC_SIF"
(11) RSEQC_REFS: "$RSEQC_REFS"
Other variables:
-q "$MIN_MAP_QUAL" 		minimum mapping quality for read to be counted (not used by all RSeQC tools)
RSeQC version: "$RSEQC_VERSION"
"

EXPECTED_ARGS=11
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
    echo -e "Usage: "$SCRIPT_TITLE" <SEQ_TYPE> <ANALYSIS_DIR> <QC_DIR> <SAMPLE_DIR> <ALIGNMENT_DIRNAME> <BAM_IN_FILENAME> <SAMPLE_NAME> <REFSEQ_BED> <HOUSEKEEPING_BED> <RSEQC_SIF> <RSEQC_REFS>
    ${red}ERROR - expecting "$EXPECTED_ARGS" ARGS but "$#" were provided:${NC}
    "$@"
    "
    exit 1
fi

# IF statements to check/create directories
	# # Check for Sample Alignment directory
	# if [ ! -d "$SAMPLE_DIR"/Alignments ]
	# then
	# 	echo -e "${red}ERROR - folder does not exist: "$SAMPLE_DIR"/Alignments
	# 	"
	# 	exit 1
	
	# elif [ -d "$SAMPLE_DIR"/Alignments ]
	# then
	# 	echo "Alignments folder found: "$SAMPLE_DIR"/Alignments
	# 	"
	# fi

	# Check for BAM file in alignment directory
	BAM_IN_FILE=$SAMPLE_DIR/$ALIGNMENT_DIRNAME/$BAM_IN_FILENAME 
	if [ ! -f $BAM_IN_FILE ]
	then
		echo "${red}ERROR - file does not exist: $BAM_IN_FILE"
		echo
		exit 1
	elif [ -f $BAM_IN_FILE ]
	then
		echo "BAM file found: $BAM_IN_FILE"
	fi

	# Check for QC directories
	OUTPUT_DIR="$QC_DIR"/RSeQC_metrics/
	if [ ! -d "$OUTPUT_DIR" ]
	then
		echo "Creating directory: "$OUTPUT_DIR"
		"
		mkdir --parents "$OUTPUT_DIR"

	else
		echo "Output directory found: "$OUTPUT_DIR""
	fi

# BAMFILE2="$SAMPLE_DIR"/Alignments/"$SAMPLE_NAME".mapped.rgid.filtered.sorted.dups_mark.bam
OUTPUT_FILE_BASE="$OUTPUT_DIR"/"$SAMPLE_NAME"_RSeQC

cd "$OUTPUT_DIR"

echo -e "
${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"



	echo "
	Running bam_stat.py..."
	(>&2 echo "
		Error messages for bam_stat.py:")
	singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RSEQC_REFS" "$RSEQC_SIF" bam_stat.py \
		-i "$BAM_IN_FILE" \
		-q "$MIN_MAP_QUAL" | tee "$OUTPUT_FILE_BASE".bam_stat.txt

		# check output status
		if [ $? -ne 0 ]
		then
			echo -e "${red}bam_stat.py failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi


	echo "
	Running geneBody_coverage.py..."
	(>&2 echo "
		Error messages for geneBody_coverage.py:")
	singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RSEQC_REFS" "$RSEQC_SIF" geneBody_coverage.py \
		-i "$BAM_IN_FILE" \
		-r "$HOUSEKEEPING_BED" \
		-f pdf \
		-o "$OUTPUT_FILE_BASE"
	
		# check output status
		if [ $? -ne 0 ]
		then
			echo -e "${red}geneBody_coverage.py failed${NC}
			"
			exit 1
		else
			echo "Done."
			rm log.txt 	# File not needed
		fi


	echo "
	Running infer_experiment.py..."
	(>&2 echo "
		Error messages for infer_experiment.py:")
	singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RSEQC_REFS" "$RSEQC_SIF" infer_experiment.py \
		-i "$BAM_IN_FILE" \
		-r "$REFSEQ_BED" \
		-s 1000000 | tee "$OUTPUT_FILE_BASE".infer_experiment.txt

		# check output status
		if [ $? -ne 0 ]
		then
			echo -e "${red}infer_experiment.py failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi


	if [ "$SEQ_TYPE" = "PE" ]
		then	
			echo "Running inner_distance.py..."
			(>&2 echo "
				Error messages for inner_distance.py:")
			singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RSEQC_REFS" "$RSEQC_SIF" inner_distance.py \
				-i "$BAM_IN_FILE" \
				-o "$OUTPUT_FILE_BASE" \
				-r "$REFSEQ_BED" \
				-k 10000000 \
				-l -500 -u 500 -s 5 \
				-q "$MIN_MAP_QUAL" | tee "$OUTPUT_FILE_BASE".inner_distance_stats.txt
			
				# check output status
				if [ $? -ne 0 ]
				then
					echo -e "${red}inner_distance.py failed${NC}
					"
					exit 1
				else
					echo "Done."
					rm "$OUTPUT_FILE_BASE".inner_distance.txt  # This file just contains info about the pairs used
				fi
	fi


	echo "
	Running junction_annotation.py..."
	(>&2 echo "
			Error messages for junction_annotation.py:")
	singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RSEQC_REFS" "$RSEQC_SIF" junction_annotation.py \
		-i "$BAM_IN_FILE" \
		-r "$REFSEQ_BED" \
		-o "$OUTPUT_DIR"/"$SAMPLE_NAME".rseqc.junction_annotation \
		-q "$MIN_MAP_QUAL"
		
		# check output status
		if [ $? -ne 0 ]
		then
			echo -e "${red}junction_annotation.py failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi


	echo "
	Running junction_saturation.py..."
	(>&2 echo "
			Error messages for junction_saturation.py:")
	singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RSEQC_REFS" "$RSEQC_SIF" junction_saturation.py \
		-i "$BAM_IN_FILE" \
		-r "$REFSEQ_BED" \
		-o "$OUTPUT_FILE_BASE" \
		-q "$MIN_MAP_QUAL"
		
		# check output status
		if [ $? -ne 0 ]
		then
			echo -e "${red}junction_saturation.py failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi


	echo "
	Running read_distribution.py..."
	(>&2 echo "
			Error messages for read_distribution.py:")
	singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RSEQC_REFS" "$RSEQC_SIF" read_distribution.py \
		-i "$BAM_IN_FILE" \
		-r "$REFSEQ_BED" | tee "$OUTPUT_FILE_BASE".read_distribution.txt
	    # “Total Reads”: This does NOT include those QC fail,duplicate and non-primary hit reads
	    # “Total Tags”: reads spliced once will be counted as 2 tags, reads spliced twice will be counted as 3 tags, etc. And because of this, “Total Tags” >= “Total Reads”
	    # “Total Assigned Tags”: number of tags that can be unambiguously assigned the 10 groups (see below table).
	    # Tags assigned to “TSS_up_1kb” were also assigned to “TSS_up_5kb” and “TSS_up_10kb”, tags assigned to “TSS_up_5kb” were also assigned to “TSS_up_10kb”. Therefore, “Total Assigned Tags” = CDS_Exons + 5’UTR_Exons + 3’UTR_Exons + Introns + TSS_up_10kb + TES_down_10kb.
	    # When assign tags to genome features, each tag is represented by its middle point.
		
		# check output status
		if [ $? -ne 0 ]
		then
			echo -e "${red}read_distribution.py failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi


	echo "
	Running read_duplication.py..."
	(>&2 echo "
			Error messages for read_duplication.py:")
	singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RSEQC_REFS" "$RSEQC_SIF" read_duplication.py \
	-i "$BAM_IN_FILE" \
	-o "$OUTPUT_FILE_BASE" \
	-q "$MIN_MAP_QUAL"
	# Two strategies are used to determine reads duplication rate: * Sequence based: reads with identical sequence are regarded as duplicated reads. * Mapping based: reads mapped to the exactly same genomic location are regarded as duplicated reads. For splice reads, reads mapped to the same starting position and splice the same way are regarded as duplicated reads.

		# check output status
		if [ $? -ne 0 ]
		then
			echo -e "${red}read_duplication.py failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi


	echo "
	Running read_GC.py..."
	(>&2 echo "
			Error messages for read_GC.py:")
	singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RSEQC_REFS" "$RSEQC_SIF" read_GC.py \
		-i "$BAM_IN_FILE" \
		-o "$OUTPUT_FILE_BASE" \
		-q "$MIN_MAP_QUAL"
		
		# check output status
		if [ $? -ne 0 ]
		then
			echo -e "${red}read_GC.py failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi


	echo "
	Running read_NVC.py..."
	(>&2 echo "
			Error messages for read_NVC.py:")
	singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RSEQC_REFS" "$RSEQC_SIF" read_NVC.py \
		-i "$BAM_IN_FILE" \
		-o "$OUTPUT_FILE_BASE" \
		-q "$MIN_MAP_QUAL"
		
		# check output status
		if [ $? -ne 0 ]
		then
			echo -e "${red}read_NVC.py failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi

# This does not appear to work but does not report any errors and therefore causes problems for the pipeline
	# echo "
	# Running read_quality.py..."
	# (>&2 echo "
	# 		Error messages for read_quality.py:")
	# read_quality.py -i "$BAMFILE2" -o "$OUTPUT_FILE_BASE" -q "$MIN_MAP_QUAL"
		
	# 	# check output status
	# 	if [ $? -ne 0 ]
	# 	then
	# 		echo -e "${red}read_quality.py failed${NC}
	# 		"
	# 		exit 1
	# 	else
	# 		echo "Done."
	# 	fi


	echo "
	Running tin.py..."
	(>&2 echo "
			Error messages for tin.py:")
	singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RSEQC_REFS" "$RSEQC_SIF" tin.py \
	-i "$BAM_IN_FILE" \
	-r "$HOUSEKEEPING_BED"
	# This program is designed to evaluate RNA integrity at transcript level. TIN calculates a score (0 <= TIN <= 100) for each expressed transcript, however, the medTIN (i.e. meidan TIN score across all the transcripts) can also be used to measure the RNA integrity at sample level.
		
		# check output status
		if [ $? -ne 0 ]
		then
			echo -e "${red}tin.py failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi

echo -e "
${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}
"
