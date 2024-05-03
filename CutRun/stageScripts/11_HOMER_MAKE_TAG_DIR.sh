#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=11_HOMER_MAKE_TAG_DIR.sh
SCRIPT_VERSION=v0.3
# DATE: 11-23-2022
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample and is executed by <XXseq>_pipeline.sh
# This version (0.2) uses Homer singularity container.
# v0.3 050224: updating to run on Proton2

# HOMER can read in from BAM/SAM or BED format
# If your alignment is in a different format, it is recommended that you convert it into a BED file format since this is a fairly standard format with minimal chance of misinterpretation:
# Column1: chromosome
# Column2: start position
# Column3: end position
# Column4: Name (or strand +/-)
# Column5: Ignored**
# Column6: Strand +/-
# **Unfortunately, BED files are used in different ways by different groups.  By default, HOMER ignores the 5th column of BED files.  In some cases the 5th column is used to encode the number of reads at that position.  If this is the case, use the "-force5th".

# HOMER *.tags.tsv files
# These are files used to store sequencing data in HOMER tag directories.  They are tab-delimited text files that are sorted to allow for relatively quick access and processing. 
# 1. blank (can be used for a name)
# 2. chromsome
# 3. position (1-indexed)
# 4. strand (0 or 1, +/- not allowed here)
# 5. Number of reads (can be fractional)
# 6. length of the read (optional)

# HOMER can use paired-end reads.  If you use paired-end reads for ChIP-Seq or RNA-Seq, HOMER will treat each half of the read separately (and count each as 0.5 reads), which works well for a number of applications.  If you are using stranded paired-end reads, make sure to specify "-sspe" so that HOMER will correctly interpret the intended strand for the 2nd read in the mate-pair.
# MG note: this appears to be using only read length and not template length for PE reads. This creates either gaps or overlaps in coverage.
# Should be able to overcome this limitation by writing a BED file using the TLEN (field 9) from BAM to adjust length of Read 1 from each pair and treating as single-end data.

# Something like this (for fwd strand mate 1; -f81 -F4 for rev strand mate 1):
# samtools view -f65 -F20 HCT116_WT_DMSO_3.mapped.rgid.filtered.sorted.dups_mark.bam | grep YT:Z:CP | cut -f3,4,9 | sed -e 's/-//g' | awk 'BEGIN {OFS = "\t" } NR == 1 { $4 } NR >= 1 { $4 = $2 + $3 } 1' | awk 'BEGIN {FS=OFS="\t" } $2 { $2=$2-1 } $3 { $3=$4 } $4 { $4="." } { $5="0" } { $6="+" } {print}' | head
### Subtracted 1 from field 2 to convert to 0-based, half-open (sam is 1-based, closed) = subtract 1 from field 2 and nothing from field 2? (BED format end is non-inclusive ie half-open)
# Chr   Start   End  Name Score Strand
# chr1	785124	785358	.	0	+
# chr1	785124	785358	.	0	+
# chr1	6857064	6857284	.	0	+
# chr1	7195006	7195141	.	0	+
# chr1	8072362	8072507	.	0	+
# chr1	8921357	8921451	.	0	+
# chr1	9633517	9633842	.	0	+
# chr1	9633517	9633842	.	0	+
# chr1	9633517	9633842	.	0	+
# chr1	10002681	10002810	.	0	+

## NOTE this will ignore duplicates marking for now ## (add 1024 to F to exclude marked dups)


# variables from command line via <XXseq>_pipeline.sh:
SAMPLE_DIR=${1}
HOMER_DIRNAME=${2}
BAM_IN_FILENAME=${3}
SAMPLE_NAME=${4}
STRAND_TYPE=${5}
SEQ_TYPE=${6}
REF=${7}
FRAG_LENGTH=${8}
TAGS_PER_POSITION=${9}
HOMER_SIF=${10}
HOMER_DATA=${11}
# other variables:
# HOMERPLOTS=$SHARED/espipe/Aux/scripts/RNAseq/HOMER/homerPlots.R 	# This is not currently accessible
SAMTOOLS_VERSION="$(singularity run "$HOMER_SIF" samtools --version  | head -n1)"
BED_DIRNAME="$SAMPLE_DIR"/Alignments/BED
BED_FILE="$BED_DIRNAME"/"$SAMPLE_NAME".bed
HOMER_VERSION="$(singularity run "$HOMER_SIF" perl /opt/homer/configureHomer.pl -list 2>&1 | grep "HOMER" | grep "v" | cut -f3,3)"
#HOMERPLOTS=/gpfs/jespinosa/shared/Matt/PipeLinesScripts/misc/homerPlots.R

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

# Print versions and ARGS
echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
Arguments for "$SCRIPT_TITLE":
(1) SAMPLE_DIR: "$SAMPLE_DIR"
(2) HOMER_DIRNAME: "$HOMER_DIRNAME"
(3) BAM_IN_FILENAME: "$BAM_IN_FILENAME"
(4) SAMPLE_NAME: "$SAMPLE_NAME"
(5) STRAND_TYPE: "$STRAND_TYPE"
(6) SEQ_TYPE: "$SEQ_TYPE"
(7) REF: "$REF"
(8) FRAG_LENGTH: "$FRAG_LENGTH"
(9) TAGS_PER_POSITION: "$TAGS_PER_POSITION"
(10) HOMER_SIF: "$HOMER_SIF"
(11) HOMER_DATA: "$HOMER_DATA"
Other variables:
BED_DIRNAME="$SAMPLE_DIR"/Alignments/BED
BED_FILE="$BED_DIRNAME"/"$SAMPLE_NAME".bed
NOT USED HOMERPLOTS: "$HOMERPLOTS"
Samtools version: "$SAMTOOLS_VERSION"
HOMER version: "$HOMER_VERSION"
HOMER makeTagDirectory options:
-format sam
-flip (only if STRAND_TYPE = "strand-specific-rev")
-checkGC
"

EXPECTED_ARGS=11
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo "Usage: "$SCRIPT_TITLE" <SAMPLE_DIR> <HOMER_DIRNAME> <BAM_IN_FILENAME> <SAMPLE_NAME> <STRAND_TYPE (strand-specific-fwd/strand-specific-rev/unstranded)> <SEQ_TYPE (PE/SE)> <REF> <FRAG_LENGTH> <TAGS_PER> <HOMER_SIF> <HOMER_DATA>
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

	# Check for Sample alignment BAM file
	if [ ! -f "$BAM_IN_FILENAME" ]
	then
		echo "${red}ERROR - file does not exist: "$SAMPLE_DIR"/"$ALIGNMENT_DIRNAME"/"$BAM_IN_FILENAME"${NC}
		"
		exit 1
	fi

	# # Check for existing or create new BED directory
	# if [ -d "$BED_DIRNAME" ]
	# then
	# 	echo "Removing existing BED directory: "$BED_DIRNAME""	# Should really replace to use existing BED files to save time
	# 	rm -rf "$BED_DIRNAME"									# Would require switch for BED file creation step below
	# 	echo "Making new BED directory: "$BED_DIRNAME""
	# 	mkdir --parents "$BED_DIRNAME"
	# else
	# 	echo "Making BED directory: "$BED_DIRNAME""
	# 	mkdir --parents "$BED_DIRNAME"
	# fi

	# Check for existing or create new HOMER Tag directory
	if [ -d "$HOMER_DIRNAME"/"$SAMPLE_NAME"_TagDirectory_frag-"$FRAG_LENGTH"_tbp-"$TAGS_PER_POSITION" ]
	then
		echo "Removing existing tag directory: "$HOMER_DIRNAME"/"$SAMPLE_NAME"_TagDirectory_frag-"$FRAG_LENGTH"_tbp-"$TAGS_PER_POSITION""
		rm -rf "$HOMER_DIRNAME"/"$SAMPLE_NAME"_TagDirectory_frag-"$FRAG_LENGTH"_tbp-"$TAGS_PER_POSITION"
		echo "Making new tag directory: "$HOMER_DIRNAME"/"$SAMPLE_NAME"_TagDirectory_frag-"$FRAG_LENGTH"_tbp-"$TAGS_PER_POSITION""
		mkdir --parents "$HOMER_DIRNAME"/"$SAMPLE_NAME"_TagDirectory_frag-"$FRAG_LENGTH"_tbp-"$TAGS_PER_POSITION"
	else
		echo "Making tag directory: "$HOMER_DIRNAME"/"$SAMPLE_NAME"_TagDirectory_frag-"$FRAG_LENGTH"_tbp-"$TAGS_PER_POSITION""
		mkdir --parents "$HOMER_DIRNAME"/"$SAMPLE_NAME"_TagDirectory_frag-"$FRAG_LENGTH"_tbp-"$TAGS_PER_POSITION"
	fi

# IF statements for setting makeTagDirectory options
	# Set or ignore -fragLength option
	if [ "$FRAG_LENGTH" != "none" ]
	then
		echo "Setting -fragLength "$FRAG_LENGTH""
		FRAG_LENGTH_OPTS="-fragLength"
		LENGTH="$FRAG_LENGTH"
	else
		echo "Not setting -fraglength; HOMER will use default"
		FRAG_LENGTH_OPTS=
		LENGTH=
	fi

	# Set or ignore -flip option
	if [ $STRAND_TYPE = "strand-specific-rev" ]
	then
		echo "Setting -flip for strand-specific-rev data"
		FLIP="-flip"
	else
		echo "Not setting -flip"
		FLIP=
	fi

	# Set or ignore -sspe option for stranded PE vs unstranded/SE data
	if [ $SEQ_TYPE = "PE" ]
	then
		if [ $STRAND_TYPE != "unstranded" ]
		then
			echo "Setting -sspe for strand-specific PE data"
			SSPE="-sspe"
		fi
	else
		echo "Not setting -sspe for unstranded or SE data"
		SSPE=
	fi

	# Set or ignore -tbp option
	if [ "$TAGS_PER_POSITION" != "all" ]
	then
		echo "Setting -tpb to "$TAGS_PER_POSITION""
		TBP="-tbp"
		TAGS="$TAGS_PER_POSITION"

	elif [ "$TAGS_PER_POSITION" = "all" ]
	then
		echo "Not setting -tbp; HOMER will use default"
		TBP=
		TAGS=
	fi

cd $SAMPLE_DIR

echo -e "${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"

## removed for SE version
# 	# Make BED file
# 	echo "Making BED file using TLEN (insert sizes)..."
# 	# read 1 fwd strand:
# 	srun samtools view -f65 -F20 "$BAM_IN_FILENAME" | grep YT:Z:CP | cut -f3,4,9 | sed -e 's/-//g' | awk 'BEGIN {OFS = "\t" } NR == 1 { $4 } NR >= 1 { $4 = $2 + $3 } 1' | awk 'BEGIN {FS=OFS="\t" } $2 { $2=$2-1 } $3 { $3=$4 } $4 { $4="." } { $5="0" } { $6="+" } {print}' > "$BED_FILE"
# 	# read 1 rev strand
# 	srun samtools view -f81 -F4 "$BAM_IN_FILENAME" | grep YT:Z:CP | cut -f3,4,9 | sed -e 's/-//g' | awk 'BEGIN {OFS = "\t" } NR == 1 { $4 } NR >= 1 { $4 = $2 + $3 } 1' | awk 'BEGIN {FS=OFS="\t" } $2 { $2=$2-1 } $3 { $3=$4 } $4 { $4="." } { $5="0" } { $6="-" } {print}' >> "$BED_FILE"
# 	echo -e "Done.\nCompressing BED file with gzip..."
# 	gzip "$BED_FILE"

# 	# check output status
# 	if [ $? -ne 0 ]
# 	then
# 		echo -e "${red}BED file generation failed${NC}
# 		"
# 		exit 1
# 	fi

# echo "Done making BED file for "$SAMPLE_NAME"."

	# Make HOMER tag directory:
	echo "Making HOMER Tag Directory for "$SAMPLE_NAME"..."
	srun singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$HOMER_DATA":/opt/homer "$HOMER_SIF" makeTagDirectory \
	                "$HOMER_DIRNAME"/"$SAMPLE_NAME"_TagDirectory_frag-"$FRAG_LENGTH"_tbp-"$TAGS_PER_POSITION" \
			-format sam
			"$FRAG_LENGTH_OPTS" "$LENGTH" \
			"$FLIP" \
			"$SSPE" \
			-genome "$REF" \
			-checkGC \
			"$TBP" "$TAGS"\
			"$BAM_IN_FILENAME" 2> "$HOMER_DIRNAME"/"$SAMPLE_NAME"_TagDirectory_frag-"$FRAG_LENGTH"_tbp-"$TAGS_PER_POSITION"/makeTagDirectory.out
	
	# check output status
	if [ $? -ne 0 ]
	then
		echo "${red}HOMER make TagDirectory failed${NC}
		"
		exit 1
	fi

echo "Done making HOMER Tag Directory for "$SAMPLE_NAME"."
# echo "Making HOMER plots ..."
	# Rscript $HOMERPLOTS <input_folder> <output_folder> <output_prefix>
	# Rscript "$HOMERPLOTS" "$HOMER_DIRNAME"/"$SAMPLE_NAME"_TagDirectory "$HOMER_DIRNAME"/"$SAMPLE_NAME"_TagDirectory_frag-"$FRAG_LENGTH"_tbp-"$TAGS_PER_POSITION" "$SAMPLE_NAME"

	# # check output status -does not work for Rscript
	# if [ $? -ne 0 ]
	# then
	# 	echo "HOMER plots failed
	# 	"
	# 	exit 1
	# fi

# echo "Done.
# "
echo -e "${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}
"

