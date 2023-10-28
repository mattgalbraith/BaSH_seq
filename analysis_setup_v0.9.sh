#!/bin/bash

SCRIPT_TITLE=analysis_setup.sh
SCRIPT_VERSION=0.9
# DATE: 06-30-2021
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run from Project/analysis_date/ and will generate:
# /SAMPLE/Raw directories with symbolic links to fastq files
# and
# /scripts/submitAll.sh with one command per sample in the form:
# nohup sh PATH/TO/PIPELINE_SCRIPT.sh SAMPLE_NAME ANALYSIS_DIR/ USER_NAME START_AT_STAGE END_AT_STAGE &> PIPELINE_SCRIPT.sh.SAMPLE_NAME.log &
# This script uses read1 / read2 in the third column of sample_locations.txt to deal with paired-end vs. single-end data
#
# Version 0.5: Generates new form of submit scripts that run as sbatch rather than shell jobs (should improve robustness as SLURM rather than login node will be manging pipeline script).
# MUST be used with corresponding version of pipeline script with stages called via sbatch with -W and --wrap options.
# Version 0.6: Adds SUBMIT_LOG arg passed to pipeline script for getting JOB_ID and error checking
# Version 0.7: added SBATCH --account for pipeline wrapper script; added SBATCH -p longrun and SBATCH --time=100:00:00 to prevent pipeline script from timing out (especially while waiting for resources)
# Version 0.8: Intitial setup for running on HDC Eureka
# Version 0.9 102623: updating to run on Proton2
#
# The submitAll.sh script should now be run from the appropriate Project/analysis_date/scripts directory as follows:
# sbatch -o submitAll.sh.log submitAll.sh


##### ARGUMENTS PASSED FROM COMMAND LINE - DO NOT EDIT THIS SECTION #####

# variables from command line:
PIPELINE_SCRIPT=${1}            # Path to pipeline scripts directory e.g TTseq_pipeline_v0.5.sh
START_AT_STAGE=${2}
END_AT_STAGE=${3}
# Other variables:
# THIS_USER_ACCOUNT="$(whoami)"
#if adding any new ARGS, remember to add below to output for submitAll.sh

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color

EXPECTED_ARGS=3
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo -ne "${red}ERROR - expecting "$EXPECTED_ARGS" ARGS but "$#" provided:${NC}
        "$@""
        echo "Usage: sh <PATH/TO/"$SCRIPT_TITLE"> <PATH/TO/PIPELINE_SCRIPT> <START_AT_STAGE> <END_AT_STAGE>"
        echo -e "EXAMPLE: sh /gpfs/jespinosa/shared/Matt/PipeLinesScripts/analysis_setup.sh /gpfs/jespinosa/shared/Matt/PipeLinesScripts/TTseq/TTseq_pipeline_v0.1.sh 0 8"
        exit 1
fi

# Check for sample_locations.txt
if [ ! -e sample_locations.txt ]
        then
                echo -e "${red}ERROR - sample_locations.txt not found
                "
                exit 1
fi

# Check for pipeline script
if [  ! -e "$PIPELINE_SCRIPT" ]
        then
                echo -e "${red}ERROR - Script does not exist:${NC} "$PIPELINE_SCRIPT"
                "
                exit 1
fi

# # Set THIS_USER_ACCOUNT
# if [ "$ACCOUNT" = ESP ]
#         then
#                 THIS_USER_ACCOUNT=espinosalab-"$(whoami)"
# elif [ "$ACCOUNT" = HTP ]
#         then
#                 THIS_USER_ACCOUNT=htp-"$(whoami)"
# else
#         echo -e "${red}ERROR - Incorrect account given:${NC} use either ESP or HTP
#                 "
#                 exit 1
# fi

## Set up SAMPLE directories and generate symbolic links to original raw FASTQ files
# The symbolic links allow us to use informative names without renaming or duplicating original files
# This currently requires adding a third field so sample_locations.txt with read1/read2 to deal with paired-end data (single read data will just be labelled as read1)
# New in v0.2: Will not over-write existing Sample directories or symbolic links so that analysis_setup.sh can easily be re-run to create submitAll.sh scripts for different subsets of stages.

# Create SAMPLE directory:
for i in $(cat sample_locations.txt | awk '{print $1}' | uniq)
do
        SAMPLE_NAME=$(echo "$i")
        if [ -d Sample_"$SAMPLE_NAME"/Raw ]
                then
                        echo "Raw data directory for "$SAMPLE_NAME" already present - no action needed"
                else
                        mkdir -p Sample_"$SAMPLE_NAME"/Raw
                        echo "Making directory: $(pwd)/Sample_"$SAMPLE_NAME"/Raw"
        fi
done

# Create symbolic links for read 1:
for i in $(cat sample_locations.txt | grep -w read1 | tr "\t" ":")
do
        SAMPLE_NAME=$(echo "$i" | cut -d ":" -f1)
        FASTQ_FILE=$(echo "$i" | cut -d ":" -f2)
        # Check for FASTQ_FILE
        if [ -r "$FASTQ_FILE" ]
        then
                if [ -L Sample_"$SAMPLE_NAME"/Raw/"$SAMPLE_NAME"_R1.fastq.gz ]
                then
                        echo "Symbolic link to "$SAMPLE_NAME" read 1 FASTQ file already present - no action needed"
                else
                        ln -s "$FASTQ_FILE" Sample_"$SAMPLE_NAME"/Raw/"$SAMPLE_NAME"_R1.fastq.gz
                fi
        else
                echo -e "${red}ERROR - read 1 FASTQ file does not exist or is not readable:${NC} "$FASTQ_FILE"
                "
                exit 1
        fi
done

# Create symbolic links for read 2:
for i in $(cat sample_locations.txt | grep -w read2 | tr "\t" ":")
do
        SAMPLE_NAME=$(echo "$i" | cut -d ":" -f1)
        FASTQ_FILE=$(echo "$i" | cut -d ":" -f2)
        # Check for FASTQ_FILE
        if [ -r "$FASTQ_FILE" ]
        then
                if [ -L Sample_"$SAMPLE_NAME"/Raw/"$SAMPLE_NAME"_R2.fastq.gz ]          # this does not check if linking path is same
                then
                        echo "Symbolic link to "$SAMPLE_NAME" read 2 FASTQ file already present - no action needed"
                else
                        ln -s "$FASTQ_FILE" Sample_"$SAMPLE_NAME"/Raw/"$SAMPLE_NAME"_R2.fastq.gz
                fi
        else
                echo -e "${red}ERROR - read 2 FASTQ file does not exist or is not readable:${NC} "$FASTQ_FILE"
                "
                exit 1
        fi
done

# Check for scripts directory
if [ ! -d scripts ]
then
        echo "Making directory: $(pwd)/scripts
        "
        mkdir -p $(pwd)/scripts
fi

echo "Writing commands to $(pwd)/scripts/submitAll_"$START_AT_STAGE"-"$END_AT_STAGE".sh"

# write parameters to submitAll.sh
echo "#!/bin/sh

#SBATCH -J SUBMIT_$(basename ${PIPELINE_SCRIPT%.sh})
#SBATCH -p defq
#SBATCH -o ./submitAll_"$START_AT_STAGE"-"$END_AT_STAGE".sh.log
#SBATCH -D ./
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

SCRIPT_TITLE="submitAll_"$START_AT_STAGE"-"$END_AT_STAGE".sh_for_$(basename $PWD)"
SCRIPT_VERSION="$SCRIPT_VERSION"
# DATE: $(date)
# AUTHOR: ${SCRIPT_TITLE%.sh}_v${SCRIPT_VERSION}.sh
# SUMMARY:
# This is an auto-generated script that contains commands to run the specified pipeline for each sample.
# 
# Sequencing details and analysis settings should be entered/changed in top section of <Pipeline_Type>.sh before running this script
# 
# This script should be run from the appropriate Project/analysis_date/scripts directory as follows:
# sbatch submitAll_"$START_AT_STAGE"-"$END_AT_STAGE".sh
###### ^^ COPY AND PASTE COMMAND ABOVE TO RUN ^^ #######
########################################################

# Arguments set from command line:
# Analysis Setup script: "$0"
# PIPELINE_SCRIPT: "$1"
# START_AT_STAGE: "$2"
# END_AT_STAGE: "$3"

# Arguments passed to "$PIPELINE_SCRIPT":
# THIS_SAMPLE_NAME: $(cat sample_locations.txt | awk '{print $1}' | uniq | tr "\n" ";" )
# THIS_ANALYSIS_DIR: $(pwd)/
# RAW_DIR: $(for i in $(cat sample_locations.txt | awk '{print $2}'); do dirname $i; done | uniq | tr "\n" ";")
# START_AT_STAGE: "$START_AT_STAGE"
# END_AT_STAGE: "$END_AT_STAGE"
# SUBMIT_LOG: "scripts/$(basename "$PIPELINE_SCRIPT").SAMPLE.Stage"$START_AT_STAGE"-"$END_AT_STAGE".log" 

" > scripts/submitAll_"$START_AT_STAGE"-"$END_AT_STAGE".sh

# write commands to submitAll.sh
for i in $(cat sample_locations.txt | awk '{print $1}' | uniq) # the uniq step ensures a single command per SAMPLE for PE data (rather than 1 per fastq file as listed in sample_locations.txt) (v0.2: removed the sort to preserve sensible order)
do 
        SUBMIT_LOG="$(pwd)/scripts/$(basename "$PIPELINE_SCRIPT")."$i".Stage"$START_AT_STAGE"-"$END_AT_STAGE".log"
        # NEW: NEED TO GET RAW DIR FOR MOUNTING TO SOME CONTAINERS, eg FASTQC
        RAW_DIR=$(dirname $(grep "$i" sample_locations.txt | awk '{print $2}' | head -n1))
        echo "(nohup sh "$PIPELINE_SCRIPT" "$i" "$PWD" "$RAW_DIR" "$START_AT_STAGE" "$END_AT_STAGE" "$SUBMIT_LOG" &> "$SUBMIT_LOG") &"
done >> scripts/submitAll_"$START_AT_STAGE"-"$END_AT_STAGE".sh

# need to execute commands above as subshells: (CMD) &
# then use "wait" which waits for child processes to complete
# then gather all std outs and errs to one log file:
echo "
wait

# Intialize log file:
echo  > $(pwd)/scripts/$(basename "$PIPELINE_SCRIPT").ALL.Stage"$START_AT_STAGE"-"$END_AT_STAGE".log

# Collect individual std out and err files into one log file:" >> scripts/submitAll_"$START_AT_STAGE"-"$END_AT_STAGE".sh

for i in $(cat sample_locations.txt | awk '{print $1}' | uniq)
do
        SUBMIT_LOG="$(pwd)/scripts/$(basename "$PIPELINE_SCRIPT")."$i".Stage"$START_AT_STAGE"-"$END_AT_STAGE".log"
        echo "cat "$SUBMIT_LOG" >> $(pwd)/scripts/$(basename "$PIPELINE_SCRIPT").ALL.Stage"$START_AT_STAGE"-"$END_AT_STAGE".log"
        echo "rm "$SUBMIT_LOG""
done >> scripts/submitAll_"$START_AT_STAGE"-"$END_AT_STAGE".sh

echo "Done."

exit