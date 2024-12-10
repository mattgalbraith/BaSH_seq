#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=1_FASTQ_MCF_BB.sh
SCRIPT_VERSION=0.6
# DATE: 06-30-2021
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample (incl. both reads for PE) and is executed by <XXseq>_pipline.sh
# This version (BB) uses bbduk from bbtools suite instead of fastx_trimmer - should be faster.
#
# v0.5 modified call to bbduk.sh so that it relies on PATH rather than pointing to /Matt/Tools/...
# v0.6 103023: updating to use BBTOOLS and EAUTILS Singularity/Apptainers


# variables from command line via <XXseq>_pipeline.sh:
SEQ_TYPE=${1}
READ_LENGTH=${2}
SAMPLE_DIR=${3}
SAMPLE_NAME=${4}
FASTQR1_FILE=${5}
FASTQR2_FILE=${6}
QC_DIR_NAME=${7}
MIN_QUAL=${8}
MIN_SEQ_LENGTH=${9}
MIN_ADAPTER=${10}
MIN_PERC_OCCUR=${11}
MAX_PERC_DIFF=${12}
CONTAMINANTS_FASTA=${13}
BBTOOLS_SIF=${14}
EAUTILS_SIF=${15}
THIS_ANALYSIS_DIR=${16}
RAW_DIR=${17}
# other variables:
DISCARDED_READS_DIRNAME="fastq-mcf_Discarded"
TRIMMING_SUMMARY_FILENAME="fastq-mcf_trimming_summary.txt"
FASTQ_MCF_VERSION=`singularity run "$EAUTILS_SIF" bash -c 'fastq-mcf -h | grep Version'`
BBDUK_VERSION=`singularity run "$BBTOOLS_SIF" bash -c 'echo atcg | bbduk.sh in=stdin.fa out=/dev/null 2>&1 | grep Version'`
# 2 more below:
# DISCARDED_READS_DIR="$SAMPLE_DIR"/processed/"$DISCARDED_READS_DIRNAME"
# TRIMMED_READS_DIR="$SAMPLE_DIR"/processed

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for "$SCRIPT_TITLE":
(1) SEQ_TYPE: "$SEQ_TYPE"
(2) READ_LENGTH: "$READ_LENGTH"
(3) SAMPLE_DIR: "$SAMPLE_DIR"
(4) SAMPLE_NAME: "$SAMPLE_NAME"
(5) FASTQR1_FILE: "$FASTQR1_FILE"
(6) FASTQR2_FILE: "$FASTQR2_FILE"
(7) QC_DIR_NAME: "$QC_DIR_NAME"
(8) MIN_QUAL: "$MIN_QUAL"
(9) MIN_SEQ_LENGTH: "$MIN_SEQ_LENGTH" (-l option)
(10) MIN_ADAPTER: "$MIN_ADAPTER" (-s option)
(11) MIN_PERC_OCCUR "$MIN_PERC_OCCUR" (-t option)
(12) MAX_PERC_DIFF "$MAX_PERC_DIFF" (-p option)
(13) CONTAMINANTS_FASTA: "$CONTAMINANTS_FASTA"
(14) BBTOOLS_SIF: "$BBTOOLS_SIF"
(15) EAUTILS_SIF: "$EAUTILS_SIF"
(16) THIS_ANALYSIS_DIR: "$THIS_ANALYSIS_DIR"
(17) RAW_DIR: "$RAW_DIR"
bbduk.sh version: "$BBDUK_VERSION"
bbduk.sh options:
in="$FASTQR1_FILE"
in="$FASTQR2_FILE"
out=/tmp/fastq_temp/"$(basename "$FASTQR1_FILE")"
out=/tmp/fastq_temp/"$(basename "$FASTQR2_FILE")"
int=f 	Interleaved PE = false
ftr=$[${READ_LENGTH}-1] 	Force trim right = last base to keep (0-based, exclusive)
-Xmx 64G
t=8
fastq-mcf version: "$FASTQ_MCF_VERSION"
fastq-mcf options:
-S 	(Save all discarded reads to '.skip' files)
-l "$MIN_SEQ_LENGTH" 	(Minimum remaining sequence length; default = 19)
-q "$MIN_QUAL" 	(quality threshold causing base removal; default = 10)
-s "$MIN_ADAPTER" 	(Log2 scale for adapter minimum-length-match; default=2.2 ie ~4.5)
-t "$MIN_PERC_OCCUR" 	(% occurance threshold before adapter clipping; default=0.25)
-p "$MAX_PERC_DIFF" 	(Maximum adapter difference percentage; default=10)
-o 	(output file + stats go to stdout)
"

EXPECTED_ARGS=17
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo -e "Usage: "$SCRIPT_TITLE" <SEQ_TYPE> <READ_LENGTH> <SAMPLE_DIR> <SAMPLE_NAME> <FASTQR1_FILE> <FASTQR2_FILE> <QC_DIR_NAME> <MIN_QUAL> <MIN_SEQ_LENGTH> <MIN_ADAPTER> <MIN_PERC_OCCUR> <MAX_PERC_DIFF> <CONTAMINANTS_FASTA> <BBTOOLS_SIF> <EAUTILS_SIF> <THIS_ANALYSIS_DIR> <RAW_DIR>
        ${red}ERROR - expecting "$EXPECTED_ARGS" but "$#" were provided:${NC}
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

	# check for read 1 FASTQ file
	if [ ! -f "$FASTQR1_FILE" ]
	then
		echo -e "${red}ERROR - file not found:${NC} "$FASTQR1_FILE""
		exit 1
	fi

	# check if seq type is paired-end and read 2 FASTQ file
	if [ "$SEQ_TYPE" = "PE" ]
	then
		if [ ! -f "$FASTQR2_FILE" ]
		then
			echo "${red}ERROR - file not found:${NC} "$FASTQR2_FILE""
			exit 1
		fi
	fi

	# check for contaminants file
	if [ ! -f "$CONTAMINANTS_FASTA" ]
	then
		echo "${red}ERROR - file not found:${NC} "$CONTAMINANTS_FASTA""
		exit 1
	fi

	# Set output directories
	TRIMMED_READS_DIR="$SAMPLE_DIR"/Processed
	if [ ! -d "$TRIMMED_READS_DIR" ]
	then
		echo "Making directory: "$TRIMMED_READS_DIR""
		mkdir --parents "$TRIMMED_READS_DIR"
	fi

	# Make temp directory for n+1 trimmed files
	# FASTQ_TMPDIR="$TMPDIR"/fastq_temp/"$SAMPLE_NAME"
	FASTQ_TMPDIR=/tmp/fastq_temp/"$SAMPLE_NAME" # v0.6: updated for Proton2
	#
	if [ -d "$FASTQ_TMPDIR" ]
	then
		rm -R "$FASTQ_TMPDIR" 	# Better to overwrite in case error exit leaves previous version (bbduk default prevents overwrite)
	fi
	mkdir --parents "$FASTQ_TMPDIR"

	DISCARDED_READS_DIR="$SAMPLE_DIR"/Processed/"$DISCARDED_READS_DIRNAME"
	if [ ! -d "$DISCARDED_READS_DIR" ]
	then
		echo "Making directory: "$DISCARDED_READS_DIR""
		mkdir --parents "$DISCARDED_READS_DIR"
	fi

	# check for QC/fastq-mcf directory
	if [ ! -d "$QC_DIR_NAME"/"fastq-mcf"/"$SAMPLE_NAME" ]
	then
		echo "Making directory: "$QC_DIR_NAME"/"fastq-mcf"/"$SAMPLE_NAME""
		mkdir -p "$QC_DIR_NAME"/fastq-mcf/"$SAMPLE_NAME"
	fi


cd "$SAMPLE_DIR"

echo -e "
${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]${NC}"

## CHECK FOR N+1 READ AND REMOVE IF PRESENT USING BBDUK ##
	echo "Checking read lengths for "$(basename "$FASTQR1_FILE")"..."

	LENGTHS_R1="$(zcat "$FASTQR1_FILE" | awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t");} }' | awk -v k=1000000 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' | awk -F"\t" '{print $2}' | awk '{ print length}' | sort -n | uniq)"

# Single-end mode:
	if [ "$SEQ_TYPE" = "SE" ]
		then
			if [ ${LENGTHS_R1} -eq ${READ_LENGTH} ]
				then
					echo -e "Read length specified as: "$READ_LENGTH"\nRead length of "$LENGTHS_R1" found - no action required."
					# leave FASTQR1_FILE variable as is

			elif [ ${LENGTHS_R1} -eq $[${READ_LENGTH}+1] ]
				then
					echo -e "Read length specified as: "$READ_LENGTH"\nRead length of "$[${READ_LENGTH}+1]" found - trimming off n+1 base using bbduk in single-end mode..."

					# Trim off n+1 base
					srun singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RAW_DIR":"$RAW_DIR" --bind "$FASTQ_TMPDIR":"$FASTQ_TMPDIR" "$BBTOOLS_SIF" bbduk.sh \
					in="$FASTQR1_FILE" \
					out="$FASTQ_TMPDIR"/"$(basename "$FASTQR1_FILE")" \
					int=f \
					ftr=$[${READ_LENGTH}-1] \
					-Xmx64g \
					t=8
								
					# check output status
					if [ $? -ne 0 ]
					then
							echo "${red}bbduk single-end mode failed for "$(basename "$FASTQR1_FILE")"${NC}"
							exit 1
					fi

					echo "Done.
					"

					# Update FASTQR1_FILE variable
					FASTQR1_FILE="$FASTQ_TMPDIR"/"$(basename "$FASTQR1_FILE")"

			else
				echo -e "Read length specified as: "$READ_LENGTH"\n${red}ERROR - other lengths found:${NC} "$(echo "$LENGTHS_R1" | tr "\ " "\t" | awk '{print $1"-"$NF}')"\nCheck input file."
				exit 1
			fi
	fi

# Paired-end mode:
	if [ "$SEQ_TYPE" = "PE" ]
		then
			echo "Checking read lengths for "$(basename "$FASTQR2_FILE")"..."
			LENGTHS_R2="$(zcat "$FASTQR2_FILE" | awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t");} }' | awk -v k=1000000 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' | awk -F"\t" '{print $2}' | awk '{ print length}' | sort -n | uniq)"

			if [ ${LENGTHS_R1} -eq ${LENGTHS_R2} ] 
				then
					if [ ${LENGTHS_R2} -eq $[${READ_LENGTH}] ]
						then
							echo "Read length specified as: "$READ_LENGTH"
							Read lengths of "$LENGTHS_R2" found - no action required."
							# leave FASTQR2_FILE variable as is
					
					elif [ ${LENGTHS_R2} -eq $[${READ_LENGTH}+1] ]
						then
							echo -e "Read length specified as: "$READ_LENGTH"\nRead lengths of "$[${READ_LENGTH}+1]" found - trimming off n+1 base using bbduk in paired-end mode..."
							# Trim off n+1 base
							srun singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RAW_DIR":"$RAW_DIR" --bind "$FASTQ_TMPDIR":"$FASTQ_TMPDIR" "$BBTOOLS_SIF" bbduk.sh \
							in1="$FASTQR1_FILE" \
							in2="$FASTQR2_FILE" \
							out1="$FASTQ_TMPDIR"/"$(basename "$FASTQR1_FILE")" \
							out2="$FASTQ_TMPDIR"/"$(basename "$FASTQR2_FILE")" \
							int=f \
							ftr=$[${READ_LENGTH}-1] \
							-Xmx64g \
							t=8
							
							# check output status
							if [ $? -ne 0 ]
							then
									echo "${red}bbduk paired-end mode failed for "$(basename "$FASTQR1_FILE")" or "$(basename "$FASTQR2_FILE")"${NC}"
									exit 1
							fi

							echo "Done.
							"

							# Update FASTQR1_FILE and FASTQR2_FILE variables
							FASTQR1_FILE="$FASTQ_TMPDIR"/"$(basename "$FASTQR1_FILE")"
							FASTQR2_FILE="$FASTQ_TMPDIR"/"$(basename "$FASTQR2_FILE")"

					else
						echo -e "Read length specified as: "$READ_LENGTH"\n${red} ERROR - other lengths found:${NC} "$(echo "$LENGTHS_R2" | tr "\ " "\t" | awk '{print $1"-"$NF}')"\nCheck input file.
						"
						exit 1
					fi

			else
				echo -e "${red}ERROR - Read lengths for R1 and R2 do not match - check input files.${NC}
				"
				exit 1
			fi
	fi

## FASTQ-MCF ##
	FASTQR1_FILE_basename=`basename $FASTQR1_FILE`
	if [ "$SEQ_TYPE" = "SE" ]
	then
		echo "Running fastq-mcf in single-end mode for "$SAMPLE_NAME"..."

		srun singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RAW_DIR":"$RAW_DIR" --bind "$FASTQ_TMPDIR":"$FASTQ_TMPDIR" --bind "$CONTAMINANTS_FASTA":"$CONTAMINANTS_FASTA" "$EAUTILS_SIF" fastq-mcf \
			-S -q "$MIN_QUAL" -l "$MIN_SEQ_LENGTH" "$CONTAMINANTS_FASTA" "$FASTQR1_FILE" \
			-o "$TRIMMED_READS_DIR"/trimmed_"$FASTQR1_FILE_basename" 1> "$QC_DIR_NAME"/fastq-mcf/"$SAMPLE_NAME"/"$TRIMMING_SUMMARY_FILENAME"
		
		# check output status
		if [ $? -ne 0 ]
		then
				echo "fastq-mcf single-end mode failed"
				exit 1
		fi

		echo "Done."

	elif [ "$SEQ_TYPE" = "PE" ]
	then
		echo "Running fastq-mcf in paired-end mode for "$SAMPLE_NAME"..."

		FASTQR2_FILE_basename=`basename $FASTQR2_FILE`
		srun singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$RAW_DIR":"$RAW_DIR" --bind "$FASTQ_TMPDIR":"$FASTQ_TMPDIR" --bind "$CONTAMINANTS_FASTA":"$CONTAMINANTS_FASTA" "$EAUTILS_SIF" fastq-mcf \
			-S -q "$MIN_QUAL" -l "$MIN_SEQ_LENGTH" -s "$MIN_ADAPTER" -t "$MIN_PERC_OCCUR" -p "$MAX_PERC_DIFF"  "$CONTAMINANTS_FASTA" "$FASTQR1_FILE" "$FASTQR2_FILE" \
			-o "$TRIMMED_READS_DIR"/trimmed_"$FASTQR1_FILE_basename" -o "$TRIMMED_READS_DIR"/trimmed_"$FASTQR2_FILE_basename" 1> "$QC_DIR_NAME"/fastq-mcf/"$SAMPLE_NAME"/"$TRIMMING_SUMMARY_FILENAME"
		
		# check output status
		if [ $? -ne 0 ]
		then
				echo "fastq-mcf paired-end mode failed"
				exit 1
		fi

		echo "Done."

	else
		echo "${red}ERROR - SEQ_TYPE parameter not recognised:${NC} "$SEQ_TYPE""
	fi

	mv $TRIMMED_READS_DIR/*.skip.gz $DISCARDED_READS_DIR/.

	# remove trimmed FASTQ files as they should longer be needed
	if [ -d "$FASTQ_TMPDIR" ]
	then
		echo "removing temporary directory: "$FASTQ_TMPDIR""
		rm -R "$FASTQ_TMPDIR"
	fi

echo -e "${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]${NC}
"

