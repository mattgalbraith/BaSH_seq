#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=9_PE_ALIGNMENT_METRICS.sh
SCRIPT_VERSION=0.4.1
# DATE: 06-21-2018
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample and is executed by <XXseq>_pipeline.sh
# REF_FILE in .fasta format must be accompanied by .dict and .fai files created by CreateSequenceDictionary and samtools faidx, respectively; All 3 can be gzipped.
#
# Version 0.4 110923: updating to use PICARD Singularity/Apptainer on Proton2
# Version 0.4.1 240424: removing REF_FILE argument


# variables from command line via <XXseq>_pipline.sh:
ANALYSIS_DIR=${1}
QC_DIR=${2}
SAMPLE_DIR=${3}
SAMPLE_NAME=${4}
PICARD_MEM=${5}
#REF_FILE=${6} 					# required for some Picard tools; No longer in use
PICARD_SIF=${6}
SAMTOOLS_SIF=${7}
# other variables:
# SAMTOOLS_VERSION="$(samtools --version 2>&1)"
SAMTOOLS_VERSION=`singularity run "$SAMTOOLS_SIF" samtools --version | head -n1`
# CASM_VERSION="$(java -Xmx64G -jar $PICARD/picard CollectAlignmentSummaryMetrics --version 2>&1)"
CASM_VERSION=`singularity run "$PICARD_SIF" java -jar /picard.jar CollectAlignmentSummaryMetrics --version 2>&1`


blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for "$SCRIPT_TITLE":
(1) ANALYSIS_DIR: "$ANALYSIS_DIR"
(2) QC_DIR: "$QC_DIR"
(3) SAMPLE_DIR: "$SAMPLE_DIR"
(4) SAMPLE_NAME: "$SAMPLE_NAME"
(5) PICARD_MEM: "$PICARD_MEM"
(6) PICARD_SIF: "$PICARD_SIF"
(7) SAMTOOLS_SIF: "$SAMTOOLS_SIF"
Samtools version: "$SAMTOOLS_VERSION"
samtools view options:
various filters
Picard CollectAlignmentSummaryMetrics version: "$CASM_VERSION"
Picard CollectAlignmentSummaryMetrics options:
REFERENCE_SEQUENCE="$REF_FILE"
ASSUME_SORTED=true
"

EXPECTED_ARGS=8
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
    echo -e "Usage: "$SCRIPT_TITLE" <ANALYSIS_DIR> <QC_DIR> <SAMPLE_DIR> <SAMPLE_NAME> <PICARD_MEM> <PICARD_SIF> <SAMTOOLS_SIF>
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

	# Check for QC directories
	OUTPUT_DIR="$QC_DIR"/PE_metrics/
	if [ ! -d "$OUTPUT_DIR" ]
	then
		echo "Creating directory: "$OUTPUT_DIR"
		"
		mkdir --parents "$OUTPUT_DIR"

	else
		echo "Output directory found: "$OUTPUT_DIR""
	fi

	PICARD_OUTPUT_DIR="$QC_DIR"/Picard/"$SAMPLE_NAME"
	if [ ! -d "$PICARD_OUTPUT_DIR" ]
	then
		echo "Creating directory: "$PICARD_OUTPUT_DIR"
		"
		mkdir --parents "$PICARD_OUTPUT_DIR"

	else
		echo "Output directory found: "$PICARD_OUTPUT_DIR""
	fi

BAMFILE1="$SAMPLE_DIR"/Alignments/"$SAMPLE_NAME".mapped.no-rgid.bam
BAMFILE2="$SAMPLE_DIR"/Alignments/"$SAMPLE_NAME".mapped.rgid.filtered.sorted.dups_mark.bam
OUTPUT_FILE="$OUTPUT_DIR"/"$SAMPLE_NAME"_PE_alignmentMetricsSummary.txt
INSERTS_FILE="$OUTPUT_DIR"/"$SAMPLE_NAME"_TLEN_HIST.txt

cd "$ANALYSIS_DIR"

echo -e "
${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"

echo "
Paired-end alignment metrics summary for "$SAMPLE_NAME"
------------------------------------------------------------" > $OUTPUT_FILE

echo "Collecting fastq stats for "$SAMPLE_NAME"...
"
# Some basic numbers from different sources to cross-check:

echo "
FASTQC pre-trimming:" >> $OUTPUT_FILE
for i in $(ls "$QC_DIR"/FASTQC_pre_filtered/"$SAMPLE_NAME"_*_fastqc/fastqc_data.txt); do cat $i | grep -e "Filename" -e "Total Sequences" -e "Sequence length"; done | paste - - - | cut -f2- >> $OUTPUT_FILE
echo "
FASTQC post-trimming:" >> $OUTPUT_FILE
for i in $(ls "$QC_DIR"/FASTQC_post_filtered/trimmed_"$SAMPLE_NAME"_*_fastqc/fastqc_data.txt); do cat $i | grep -e "Filename" -e "Total Sequences" -e "Sequence length"; done | paste - - - | cut -f2- >> $OUTPUT_FILE

echo "Collecting alignment stats for "$SAMPLE_NAME"...
"
echo "
Unfiltered BAM: "$(basename "$BAMFILE1")"" >> $OUTPUT_FILE
echo -ne "\
Total_reads: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view "$BAMFILE1" | wc -l)"
R1: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 64 "$BAMFILE1" | wc -l)"
R2: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 128 "$BAMFILE1" | wc -l)"
Unaligned: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 5 "$BAMFILE1" | wc -l)"
Aligned: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 1 -F 260 "$BAMFILE1" | wc -l)"
Duplicates: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 1025 -F 260 "$BAMFILE1" | wc -l)"
Concordant_pairs: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 1 -F 260 "$BAMFILE1" | grep YT:Z:CP | wc -l)"
Discordant_pairs: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 1 -F 260 "$BAMFILE1" | grep YT:Z:DP | wc -l)"
alignedUnpairedMate: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 1 -F 260 "$BAMFILE1" | grep YT:Z:UP | wc -l)"
" >> $OUTPUT_FILE
# ADD MORE BOWTIE2 TAGS?

echo "
Filtered dups marked BAM: "$(basename "$BAMFILE2")"" >> $OUTPUT_FILE
echo -ne "\
Total_reads: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view "$BAMFILE2" | wc -l)"
R1: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 64 "$BAMFILE2" | wc -l)"
R2: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 128 "$BAMFILE2" | wc -l)"
Unaligned: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 5 "$BAMFILE2" | wc -l)"
Aligned: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 1 -F 260 "$BAMFILE2" | wc -l)"
Duplicates: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 1025 -F 260 "$BAMFILE2" | wc -l)"
Concordant_pairs: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 1 -F 260 "$BAMFILE2" | grep YT:Z:CP | wc -l)"
Discordant_pairs: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 1 -F 260 "$BAMFILE2" | grep YT:Z:DP | wc -l)"
alignedUnpairedMate: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 1 -F 260 "$BAMFILE2" | grep YT:Z:UP | wc -l)"
Non-dup_concordant_pairs: "$(singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 1 -F 1284 "$BAMFILE2" | grep YT:Z:UP | wc -l)"
" >> $OUTPUT_FILE
# ADD MORE BOWTIE2 TAGS?

# check output status
if [ $? -ne 0 ]
then
	echo "Collecting alignment stats failed
	"
	exit 1
fi

# TLEN histogram data for CP and DP
echo "Collecting insert sizes for "$SAMPLE_NAME"...
"

HIST_MIN=0
HIST_MAX=500
HIST_BIN_SIZE=1
TMP_CP_TLEN_FILE1="$TMPDIR"/TLEN_temp/"$SAMPLE_NAME"/CP_"$(basename "$BAMFILE1")".txt
TMP_CP_TLEN_FILE2="$TMPDIR"/TLEN_temp/"$SAMPLE_NAME"/CP_"$(basename "$BAMFILE2")".txt
TMP_DP_TLEN_FILE1="$TMPDIR"/TLEN_temp/"$SAMPLE_NAME"/DP_"$(basename "$BAMFILE1")".txt
TMP_DP_TLEN_FILE2="$TMPDIR"/TLEN_temp/"$SAMPLE_NAME"/DP_"$(basename "$BAMFILE2")".txt

mkdir --parents "$TMPDIR"/TLEN_temp/"$SAMPLE_NAME"

singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 65 -F 4 "$BAMFILE1" | grep YT:Z:CP | cut -f9 | sed -e 's/-//g' > "$TMP_CP_TLEN_FILE1"
cat "$TMP_CP_TLEN_FILE1" | R --vanilla --slave -e 'x<-scan(file="stdin", quiet=T);x.freq<-table(x);data.frame(x.freq)' | cut -d " " -f2- | tail -n +2 > "$TMPDIR"/TLEN_temp/CP_"$(basename "$BAMFILE1")"_freqTable.txt
echo "
Template length (insert size) for "$(basename "$BAMFILE1")" concordant pairs:" > "$INSERTS_FILE"
(for i in $(seq "$HIST_MIN" "$HIST_BIN_SIZE" "$HIST_MAX"); do cat "$TMPDIR"/TLEN_temp/CP_"$(basename "$BAMFILE1")"_freqTable.txt | awk -v bin=$i -v count=0 '($1==bin) {count+=$2} END {print bin,count}'; done; cat "$TMPDIR"/TLEN_temp/CP_"$(basename "$BAMFILE1")"_freqTable.txt | awk -v bin="$HIST_MAX" -v count=0 '($1>bin) {count+=$2} END {print ">"bin,count}') | tr " " "\t" >> "$INSERTS_FILE"

# check output status
if [ $? -ne 0 ]
then
	echo "${red}Collecting insert sizes failed
	"
	exit 1
fi

singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 65 -F 4 "$BAMFILE1" | grep YT:Z:DP | cut -f9 | sed -e 's/-//g' > "$TMP_DP_TLEN_FILE1"
cat "$TMP_DP_TLEN_FILE1" | R --vanilla --slave -e 'x<-scan(file="stdin", quiet=T);x.freq<-table(x);data.frame(x.freq)' | cut -d " " -f2- | tail -n +2 > "$TMPDIR"/TLEN_temp/DP_"$(basename "$BAMFILE1")"_freqTable.txt
echo "
Template length (insert size) for "$(basename "$BAMFILE1")" discordant pairs:" >> "$INSERTS_FILE"
(for i in $(seq "$HIST_MIN" "$HIST_BIN_SIZE" "$HIST_MAX"); do cat "$TMPDIR"/TLEN_temp/DP_"$(basename "$BAMFILE1")"_freqTable.txt | awk -v bin=$i -v count=0 '($1==bin) {count+=$2} END {print bin,count}'; done; cat "$TMPDIR"/TLEN_temp/DP_"$(basename "$BAMFILE1")"_freqTable.txt | awk -v bin="$HIST_MAX" -v count=0 '($1>bin) {count+=$2} END {print ">"bin,count}') | tr " " "\t" >> "$INSERTS_FILE"

# check output status
if [ $? -ne 0 ]
then
	echo "${red}Collecting insert sizes failed
	"
	exit 1
fi

singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 65 -F 4 "$BAMFILE2" | grep YT:Z:CP | cut -f9 | sed -e 's/-//g' > "$TMP_CP_TLEN_FILE2"
cat "$TMP_CP_TLEN_FILE2" | R --vanilla --slave -e 'x<-scan(file="stdin", quiet=T);x.freq<-table(x);data.frame(x.freq)' | cut -d " " -f2- | tail -n +2 > "$TMPDIR"/TLEN_temp/CP_"$(basename "$BAMFILE2")"_freqTable.txt
echo "
Template length (insert size) for "$(basename "$BAMFILE2")" concordant pairs:" >> "$INSERTS_FILE"
(for i in $(seq "$HIST_MIN" "$HIST_BIN_SIZE" "$HIST_MAX"); do cat "$TMPDIR"/TLEN_temp/CP_"$(basename "$BAMFILE2")"_freqTable.txt | awk -v bin=$i -v count=0 '($1==bin) {count+=$2} END {print bin,count}'; done; cat "$TMPDIR"/TLEN_temp/CP_"$(basename "$BAMFILE2")"_freqTable.txt | awk -v bin="$HIST_MAX" -v count=0 '($1>bin) {count+=$2} END {print ">"bin,count}') | tr " " "\t" >> "$INSERTS_FILE"

# check output status
if [ $? -ne 0 ]
then
	echo "${red}Collecting insert sizes failed
	"
	exit 1
fi

singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" "$SAMTOOLS_SIF" samtools view -f 65 -F 4 "$BAMFILE2" | grep YT:Z:DP | cut -f9 | sed -e 's/-//g' > "$TMP_DP_TLEN_FILE2"
cat "$TMP_DP_TLEN_FILE2" | R --vanilla --slave -e 'x<-scan(file="stdin", quiet=T);x.freq<-table(x);data.frame(x.freq)' | cut -d " " -f2- | tail -n +2 > "$TMPDIR"/TLEN_temp/DP_"$(basename "$BAMFILE2")"_freqTable.txt
echo "
Template length (insert size) for "$(basename "$BAMFILE2")" discordant pairs:" >> "$INSERTS_FILE"
(for i in $(seq "$HIST_MIN" "$HIST_BIN_SIZE" "$HIST_MAX"); do cat "$TMPDIR"/TLEN_temp/DP_"$(basename "$BAMFILE2")"_freqTable.txt | awk -v bin=$i -v count=0 '($1==bin) {count+=$2} END {print bin,count}'; done; cat "$TMPDIR"/TLEN_temp/DP_"$(basename "$BAMFILE2")"_freqTable.txt | awk -v bin="$HIST_MAX" -v count=0 '($1>bin) {count+=$2} END {print ">"bin,count}') | tr " " "\t" >> "$INSERTS_FILE"

# check output status
if [ $? -ne 0 ]
then
	echo "${red}Collecting insert sizes failed
	"
	exit 1
fi

rm -R "$TMPDIR"/TLEN_temp/"$SAMPLE_NAME"

# TURNING THIS OFF FOR NOW AS SLOW AND NOT FULLY WORKING
# # PICARD metrics:
# echo "Starting Picard metrics...
# "
# java -Xmx64G -jar $PICARD/picard CollectAlignmentSummaryMetrics \
# 	REFERENCE_SEQUENCE="$REF_FILE" \
# 	INPUT="$BAMFILE1" \
# 	OUTPUT="$PICARD_OUTPUT_DIR"/"$(basename "$BAMFILE1")"_AlignmentSummaryMetrics.txt \
# 	ASSUME_SORTED=true \
# 	VERBOSITY=ERROR
# 	TMP_DIR=$TMPDIR\

# # check output status
# if [ $? -ne 0 ]
# then
# 	echo "${red}PICARD CollectAlignmentSummaryMetrics failed for "$BAMFILE1"
# 	"
# 	exit 1
# fi

# java -Xmx64G -jar $PICARD/picard CollectAlignmentSummaryMetrics \
# 	REFERENCE_SEQUENCE="$REF_FILE" \
# 	INPUT="$BAMFILE2" \
# 	OUTPUT="$PICARD_OUTPUT_DIR"/"$(basename "$BAMFILE2")"_AlignmentSummaryMetrics.txt \
# 	ASSUME_SORTED=true \
# 	VERBOSITY=ERROR
# 	TMP_DIR=$TMPDIR\

# # check output status
# if [ $? -ne 0 ]
# then
# 	echo "${red}PICARD CollectAlignmentSummaryMetrics failed for "$BAMFILE2"
# 	"
# 	exit 1
# fi

# PROGRAM=null \
# PROGRAM=CollectAlignmentSummaryMetrics \
# PROGRAM=CollectInsertSizeMetrics \
# PROGRAM=CollectRnaSeqMetrics \
# PROGRAM=CollectSequencingArtifactMetrics \
# PROGRAM=CollectQualityYieldMetrics\

echo -e "${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]${NC}"

