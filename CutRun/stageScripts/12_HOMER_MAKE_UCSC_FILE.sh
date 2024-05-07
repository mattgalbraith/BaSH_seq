#!/bin/bash
# consider adding: set -e (kills script when any command returns failure code) and set -u (fails if trying to use and unset variable)

SCRIPT_TITLE=12_HOMER_MAKE_UCSC_FILE.sh
SCRIPT_VERSION=0.5
# DATE: 11-23-2022
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample and is executed by <TTseq>_pipline.sh
# This script should be run together with previous stage.
# This version (0.3) uses Homer singularity container.
# Cut & Run version - use UNSTRANDED 
# v0.4 050224: updating to run on Proton2
# v0.5 050724: adding switch for -avg argument to makeUCSCfile


# variables from command line via <XXseq>_pipline.sh:
SAMPLE_DIR=${1}
TAG_DIRNAME=${2}
SAMPLE_NAME=${3}
TRACKS_DIRNAME=${4}
#################################################################
# STRAND_TYPE=${5}
STRAND_TYPE="unstranded" # SET TO UNSTRANDED HERE FOR CUT & RUN, REPLISEQ
#################################################################
REF=${6} 					# required for making bigWigs
FRAG_LENGTH=${7}
TAGS_PER_POSITION=${8}		# used for -tbp
NORM=${9}					# used for -norm
RES=${10}					# used for -res option; default: 1; resolution in bp
THIS_ANALYSIS_DIR=${11}
HOMER_SIF=${12}
HOMER_DATA=${13}
# other variables:
HOMER_VERSION="$(singularity run "$HOMER_SIF" perl /opt/homer/configureHomer.pl -list 2>&1 | grep "HOMER" | grep "v" | cut -f3,3)"
FILESIZE=1e10 				# used for -fsize option; default: 1e10; no reduction


blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

echo "
Script name: "$SCRIPT_TITLE"
Script version: "$SCRIPT_VERSION"
HOMER version: "$HOMER_VERSION"
Arguments for "$SCRIPT_TITLE":
(1) SAMPLE DIR: "$SAMPLE_DIR"
(2) HOMER TAG DIRNAME: "$TAG_DIRNAME"
(3) SAMPLE NAME: "$SAMPLE_NAME"
(4) TRACKS DIRECTORY: "$TRACKS_DIRNAME"
(HARD CODED) STRAND_TYPE: "$STRAND_TYPE"
(6) REF: "$REF"
(7) FRAG LENGTH: "$FRAG_LENGTH"
(8) TAGS PER POSITION: "$TAGS_PER_POSITION"
(9) RESOLUTION: "$RES"
(10) NORMALIZATION TO: "$NORM"
(11) THIS_ANALYSIS_DIR: "$THIS_ANALYSIS_DIR"
(12) HOMER_SIF: "$HOMER_SIF"
(13) HOMER_DATA: "$HOMER_DATA"
HOMER version: "$HOMER_VERSION"
HOMER makeTagDirectory options:
Other options:
TAG DIRECTORY: "$TAG_DIRNAME"
FILE SIZE: "$FILESIZE"
"

EXPECTED_ARGS=13
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
    echo -e "Usage: "$SCRIPT_TITLE" <SAMPLE_DIR> <HOMER_TAG_DIRNAME> <SAMPLE_NAME> <TRACKS_DIRNAME> <STRAND_TYPE (strand-specific-fwd/strand-specific-rev/unstranded)> <REF> <FRAG_LENGTH> <TAGS_PER> <NORM> <RES> <ANALYSIS_DIR> <HOMER_SIF> <HOMER_DATA>
	    ${red}ERROR - expecting "$EXPECTED_ARGS" ARGS but "$#" were provided:${NC}
        "$@"
        "
    exit 1
fi

# IF statements to check/create directories
	# Check for Sample directory
	if [ ! -d $SAMPLE_DIR ]
	then
		echo -e "${red}ERROR - folder does not exist: "$SAMPLE_DIR"${NC}
		"
		exit 1
	fi

	# Check for Sample Tag directory
	if [ ! -d "$TAG_DIRNAME" ]
	then
		echo -e "${red}ERROR - Tag directory not found: "$TAG_DIRNAME"${NC}
		"
		exit 1
	fi

	# Check for tracks directory 						# Changed location to analysis_dir/Tracks/HOMER_frag-"$FRAG_LENGTH"_tbp-"$TAGS_PER_POSITION"
	OUTPUT_DIR="$TRACKS_DIRNAME"/HOMER_frag-"$FRAG_LENGTH"_tbp-"$TAGS_PER_POSITION"_res-"$RES"
	if [ -d "$OUTPUT_DIR" ]
		then
			echo -e "Tracks directory found: "$OUTPUT_DIR""
	elif [ ! -d "$OUTPUT_DIR" ]
		then
			echo -e "Making tracks directory: "$OUTPUT_DIR""
			mkdir --parents "$OUTPUT_DIR"
	fi

	# Remove existing tracks with same name
	if ls "$OUTPUT_DIR"/"$SAMPLE_NAME"*.bedGraph* 1> /dev/null 2>&1
	then	
		echo "Removing existing tracks for "$SAMPLE_NAME" in: "$OUTPUT_DIR""
		rm "$OUTPUT_DIR"/"$SAMPLE_NAME"*.bedGraph*
	fi

# IF statements for setting makeUCSCfile options
	# Set or ignore -fragLength option
	if [ "$FRAG_LENGTH" != "none" ]
	then
		echo "Setting -fragLength "$FRAG_LENGTH""
		FRAG_LENGTH_OPTS="-fragLength"
		LENGTH="$FRAG_LENGTH"
	else
		echo "Not using -fraglength"
		FRAG_LENGTH_OPTS=
		LENGTH=
	fi

	# Set or ignore -tbp option
	if [ "$TAGS_PER_POSITION" != "all" ]
	then
		echo "Setting -tpb to "$TAGS_PER_POSITION""
		TBP="-tbp"
		TAGS="$TAGS_PER_POSITION"

	elif [ "$TAGS_PER_POSITION" = "all" ]
	then
		echo "not setting -tbp"
		TBP=
		TAGS=
	fi

	# Set or ignore -avg option
	if [ "$RES" == 1 ]
	then
		echo "not setting -avg"
		AVG=

	elif [ "$RES" != 1 ]
	then
		echo "setting -avg"
		AVG="-avg"
	fi


cd "$SAMPLE_DIR"

echo -e "${blue}"$SCRIPT_TITLE" STARTED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]" && echo -e "${NC}"

echo "HOMER making UCSC file(s) for "$SAMPLE_NAME"..."

if [ "$STRAND_TYPE" = "unstranded" ]
	then
		echo "Making normalized unstranded bedGraph..."

		srun singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$HOMER_DATA":/opt/homerdata "$HOMER_SIF" makeUCSCfile $(echo -e "\
					 "$TAG_DIRNAME" \
					 "$FRAG_LENGTH_OPTS" "$LENGTH" \
					 -norm "$NORM" \
					 -fsize "$FILESIZE" \
					 -color 0,0,255 \
					 -name "$SAMPLE_NAME".unstranded.norm"$NORM".frag-"$FRAG_LENGTH".tbp-"$TAGS_PER_POSITION" \
					 -res "$RES" \
					 "$AVG" \
					 "$TBP" "$TAGS" \
					 -strand both \
					 ")\
					 > "$OUTPUT_DIR"/"$SAMPLE_NAME".unstranded.norm"$NORM".frag-"$FRAG_LENGTH".tbp-"$TAGS_PER_POSITION"_res-"$RES".bedGraph

		if [ $? -ne 0 ]
		then
			echo -e "${red}Making normalized unstranded bedGraph failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi

		echo "Compressing bedGraph with gzip..."
		gzip "$OUTPUT_DIR"/"$SAMPLE_NAME".unstranded.norm"$NORM".frag-"$FRAG_LENGTH".tbp-"$TAGS_PER_POSITION".bedGraph

		if [ $? -ne 0 ]
		then
			echo -e "${red}Gzip failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi
		# # # Extra steps to convert bedGraph to bigWig
		# echo "making normalized bigwig"
		# cat $WIG_DIR/$SAMPLE_NAME.normalized.bedGraph.unsorted | grep "^track" > $WIG_DIR/$SAMPLE_NAME.normalized.bedGraph
		# cat $WIG_DIR/$SAMPLE_NAME.normalized.bedGraph.unsorted | grep -v "^track" | sort -k1,1 -k2,2n >> $WIG_DIR/$SAMPLE_NAME.normalized.bedGraph
		# bedGraphToBigWig $WIG_DIR/$SAMPLE_NAME.normalized.bedGraph $REF_INDEX $WIG_DIR/$SAMPLE_NAME.normalized.bigWig
		# if [ $? -ne 0 ]; then exit 1; fi
		# gzip $WIG_DIR/$SAMPLE_NAME.normalized.bigWig
		# rm -f $WIG_DIR/$SAMPLE_NAME.normalized.bedGraph.unsorted

	elif [ "$STRAND_TYPE" != "unstranded" ]
	then
		echo "Making normalized positive-strand bedGraph..."

		srun singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$HOMER_DATA":/opt/homer "$HOMER_SIF" makeUCSCfile $(echo -e "\
					 "$TAG_DIRNAME" \
					 "$FRAG_LENGTH_OPTS" "$LENGTH" \
					 -norm "$NORM" \
					 -fsize "$FILESIZE" \
					 -color 0,0,255 \
					 -name "$SAMPLE_NAME".pos-strand.norm"$NORM".frag-"$FRAG_LENGTH".tbp-"$TAGS_PER_POSITION" \
					 -res "$RES" \
					 "$AVG" \
					 "$TBP" "$TAGS" \
			 		 -strand + \
					 ")\
					 > "$OUTPUT_DIR"/"$SAMPLE_NAME".pos-strand.norm"$NORM".frag-"$FRAG_LENGTH".tbp-"$TAGS_PER_POSITION"_res-"$RES".bedGraph

		if [ $? -ne 0 ]
		then
			echo -e "${red}Making normalized positive-strand bedGraph failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi

		echo "Compressing bedGraph with gzip..."
		gzip "$OUTPUT_DIR"/"$SAMPLE_NAME".pos-strand.norm"$NORM".frag-"$FRAG_LENGTH".tbp-"$TAGS_PER_POSITION".bedGraph

		if [ $? -ne 0 ]
		then
			echo -e "${red}Gzip failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi
	# # Extra steps to convert bedGraph to bigWig
	# echo "making normalized positive bigwig"
	# cat $WIG_DIR/$SAMPLE_NAME.pos.normalized.bedGraph.unsorted | grep "^track" > $WIG_DIR/$SAMPLE_NAME.pos.normalized.bedGraph
	# cat $WIG_DIR/$SAMPLE_NAME.pos.normalized.bedGraph.unsorted | grep -v "^track" | sort -k1,1 -k2,2n >> $WIG_DIR/$SAMPLE_NAME.pos.normalized.bedGraph
	# bedGraphToBigWig $WIG_DIR/$SAMPLE_NAME.pos.normalized.bedGraph $REF_INDEX $WIG_DIR/$SAMPLE_NAME.pos.normalized.bigWig
	# if [ $? -ne 0 ]; then exit 1; fi
	# gzip $WIG_DIR/$SAMPLE_NAME.pos.normalized.bedGraph
	# rm -f $WIG_DIR/$SAMPLE_NAME.pos.normalized.bedGraph.unsorted

		echo "Making normalized negative-strand bedGraph..."

		srun singularity run --bind "$THIS_ANALYSIS_DIR":"$THIS_ANALYSIS_DIR" --bind "$HOMER_DATA":/opt/homer "$HOMER_SIF" makeUCSCfile $(echo -e "\
					 "$TAG_DIRNAME" \
					 "$FRAG_LENGTH_OPTS" "$LENGTH" \
					 -norm "$NORM" \
					 -fsize "$FILESIZE" \
					 -color 255,0,0 \
					 -name "$SAMPLE_NAME".neg-strand.norm"$NORM".frag-"$FRAG_LENGTH".tbp-"$TAGS_PER_POSITION" \
					 -res "$RES" \
					 "$AVG" \
					 "$TBP" "$TAGS" \
					 -strand - \
					 ")\
					 > "$OUTPUT_DIR"/"$SAMPLE_NAME".neg-strand.norm"$NORM".frag-"$FRAG_LENGTH".tbp-"$TAGS_PER_POSITION"_res-"$RES".bedGraph

		if [ $? -ne 0 ]
		then
			echo -e "${red}Making normalized positive-strand bedGraph failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi

		echo "Compressing bedGraph with gzip..."
		gzip "$OUTPUT_DIR"/"$SAMPLE_NAME".neg-strand.norm"$NORM".frag-"$FRAG_LENGTH".tbp-"$TAGS_PER_POSITION"_res-"$RES".bedGraph

		if [ $? -ne 0 ]
		then
			echo -e "${red}Gzip failed${NC}
			"
			exit 1
		else
			echo "Done."
		fi
	# # Extra steps to convert bedGraph to bigWig
	# echo "making normalized negative bigwig"
	# cat $WIG_DIR/$SAMPLE_NAME.neg.normalized.bedGraph.unsorted | grep "^track" > $WIG_DIR/$SAMPLE_NAME.neg.normalized.bedGraph
	# cat $WIG_DIR/$SAMPLE_NAME.neg.normalized.bedGraph.unsorted | grep -v "^track" | sort -k1,1 -k2,2n >> $WIG_DIR/$SAMPLE_NAME.neg.normalized.bedGraph
	# bedGraphToBigWig $WIG_DIR/$SAMPLE_NAME.neg.normalized.bedGraph $REF_INDEX $WIG_DIR/$SAMPLE_NAME.neg.normalized.bigWig
	# if [ $? -ne 0 ]; then exit 1; fi
	# gzip $WIG_DIR/$SAMPLE_NAME.neg.normalized.bedGraph
	# rm -f $WIG_DIR/$SAMPLE_NAME.neg.normalized.bedGraph.unsorted

fi

echo -e "${green}"$SCRIPT_TITLE" ENDED AT: " `date` "[JOB_ID:" $SLURM_JOB_ID" NODE_NAME:" $SLURMD_NODENAME"]${NC}
"

