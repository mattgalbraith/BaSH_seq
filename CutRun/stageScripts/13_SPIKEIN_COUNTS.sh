#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=13_SPIKEIN_COUNTS.sh
SCRIPT_VERSION=v0.1
# DATE: 10-12-2023
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to count exact matches (using grep) to nucleosome spike-in barcodes in each FASTQ file.
# This script is designed to be run once for each sample and is executed by <XXseq>_pipeline.sh


# variables from command line via <XXseq>_pipeline.sh:
SEQ_TYPE=${1}
SAMPLE_NAME=${2}
FASTQR1_FILE=${3}
FASTQR2_FILE=${4}
SPIKEIN_DIRNAME=${5}
# other variables:
# none at present

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

# Print versions and ARGS
echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for 13_SPIKEIN_COUNTS.sh:
(1) SEQ_TYPE: "$SEQ_TYPE"
(2) SAMPLE_NAME: "$SAMPLE_NAME"
(3) FASTQR1_FILE: "$FASTQR1_FILE"
(4) FASTQR2_FILE: "$FASTQR2_FILE"
(5) SPIKEIN_DIRNAME: "$SPIKEIN_DIRNAME"
"

EXPECTED_ARGS=5
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo "Usage: "$SCRIPT_TITLE" <SEQ_TYPE> <SAMPLE_NAME> <FASTQR1_FILE> <FASTQR2_FILE> <SPIKEIN_DIRNAME>
        ${red}ERROR - expecting "$EXPECTED_ARGS" ARGS but "$#" were provided:${NC}
        "$@"
        "
        exit 1
fi

# IF statements to check/create directories
	# Check for Sample directory
	if [ ! -d $SAMPLE_DIR ]
	then
		echo "${red}ERROR - folder does not exist: "$SAMPLE_DIR"${NC}
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

	# Check for spike-in counts output directory
	if [ ! -d $SAMPLE_DIR/$SPIKEIN_DIRNAME ]
	then
		mkdir $SAMPLE_DIR/$SPIKEIN_DIRNAME
	fi



cd $SAMPLE_DIR

echo -e "${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"

if [ "$SEQ_TYPE" = "SE" ]
		then
			echo "Running spike-in counting in single-end mode for "$SAMPLE_NAME"...
			"

			# () are used to group commands to run in a subshell, allowing all output to be redirected to single file
			# set -o pipefail sets exit status of pipe/subshell to exit code of last command to exit non-zero (otherwise non-zero from bowtie2 will not be caught)
			# outer brackets with pipefail ensures that we do not just capture the exit code from tee
			(set -o pipefail
				(for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA
					do echo $barcode; zcat $FASTQR1_FILE | grep -c $barcode; done | paste - - > ""$SPIKEIN_DIRNAME"/"$SAMPLE_NAME"_R1_spikein_counts.txt"))
			
			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}spike-in counting failed in single-end mode"
					exit 1
			fi

	elif [ "$SEQ_TYPE" = "PE" ]
		then
			echo "Running spike-in counting in paired-end mode for "$SAMPLE_NAME"...
			"
			# () are used to group commands to run in a subshell, allowing all output to be redirected to single file
			# set -o pipefail sets exit status to exit code of pipe/subshell to last command to exit non-zero (otherwise non-zero from hisat2 will not be caught)
			# outer brackets with pipefail ensures that we do not just capture the exit code from tee
			(set -o pipefail
				(for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA
					do echo $barcode; zcat $FASTQR1_FILE | grep -c $barcode; done | paste - - > ""$SPIKEIN_DIRNAME"/"$SAMPLE_NAME"_R1_spikein_counts.txt"))

			(set -o pipefail
				(for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA
					do echo $barcode; zcat $FASTQR2_FILE | grep -c $barcode; done | paste - - > ""$SPIKEIN_DIRNAME"/"$SAMPLE_NAME"_R2_spikein_counts.txt"))
			
			# check output status
			if [ $? -ne 0 ]
			then
					echo -e "${red}spike-in counting failed in paired-end mode"
					exit 1
			fi

	else
			echo -e "${red}ERROR - SEQ_TYPE parameter not recognised:${NC} "$SEQ_TYPE"
			"
	fi

echo -e "${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}
"

