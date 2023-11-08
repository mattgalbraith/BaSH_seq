#!/bin/bash -l

#SBATCH -o ./FASTQ_MERGE_2.out
#SBATCH -e ./FASTQ_MERGE_2.err
#SBATCH -D ./
#SBATCH -J FASTQ_MERGE_2

SCRIPT_TITLE=FASTQ_MERGE_2.sh
SCRIPT_VERSION=2.1
# DATE: 11-29-2021
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be called from submit_fastq_merging.sh and will merge (concatenate) multiple fastq files per sample from multiple sequencing runs.
# Requires parameters supplied on command line:
# $1 = Desired SAMPLE_NAMES 
# $2 = read1/read2 
# $3-5 = full paths to fastq files to be merged in remaining.
# Merged file will use basename of file in 3rd column (i.e. the first fastq file).
# A sample_locations.txt file will be automatically created for merged files; this can be copied to new analysis_dir.

# directly concatenates the gzipped fastq files which when uncompressed give the same result (merged file size will be slightly larger then the longer zcat>merged.fastq>gzip method).
# see: https://stackoverflow.com/questions/8005114/fast-concatenation-of-multiple-gzip-files

echo ""$SCRIPT_TITLE" version: "$SCRIPT_VERSION""
echo "Started at: "$(date "+%Y-%m-%d-%H%M")""
echo "------"

                FILES=$(echo "$3 $4 $5")
                NAME=$(basename "$3")
                echo "Merging files for ${NAME%.fastq.gz}..."
                echo "Source files:"
                echo "$(for i in $(echo "$FILES" | tr "\t" "\n"); do echo ""$i"    Reads: $(( $(zcat $i | wc -l) / 4 ))"; done)"
                echo "Merged file: "$PWD"/${NAME%.fastq.gz}_merged.fastq.gz"
                cat $FILES > ${NAME%.fastq.gz}_merged.fastq.gz
                echo "Checking final read count for merged file:"
                echo $(( $(zcat ${NAME%.fastq.gz}_merged.fastq.gz | wc -l) / 4 ))
                echo "Writing to sample_locations.txt"
                echo -e "$1\t$PWD"/${NAME%.fastq.gz}_merged.fastq.gz"\t$2" >> ${NAME%.fastq.gz}_sample_locations.txt
                echo "Done."
        done

echo "------"
echo "FINISHED at "$(date "+%Y-%m-%d-%H%M")""