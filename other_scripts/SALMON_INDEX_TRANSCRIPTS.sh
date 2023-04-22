#!/bin/bash -l

#SBATCH -o ./SALMON_INDEX_TRANSCRIPTS.out
#SBATCH -e ./SALMON_INDEX_TRANSCRIPTS.err
#SBATCH -D ./
#SBATCH -J SALMON_INDEX_TRANSCRIPTS
#SBATCH --partition=c2s30


# DATE: 12-05-2022
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to build a Salmon index for transcript quantification of RNA-seq data, running Salmon from a Singularity/Apptainer container.
# See https://combine-lab.github.io/salmon/getting_started/
# See https://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html#salmon-quantification
# See https://bioconductor.org/packages/release/bioc/vignettes/fishpond/inst/doc/swish.html#The_Swish_method


SCRIPT_TITLE=SALMON_INDEX_TRANSCRIPTS.sh
SCRIPT_VERSION=0.1
SALMON_SIF=/home/matthew.galbraith/Tools/Singularity/salmon.sif
SALMON_VERSION="$(singularity run $SALMON_SIF salmon --version)"

cd /home/matthew.galbraith/references/GRCh38_salmon

echo "Starting "$SCRIPT_TITLE" at "$(date "+%Y-%m-%d-%H%M")""
echo "Version = $SALMON_VERSION"

singularity run $SALMON_SIF salmon index \
        -p 10 \
        --gencode \
        -t /home/matthew.galbraith/references/GRCh38/Gencode_r42/fasta/gencode.v42.transcripts.fa.gz \
        -i gencode.v42_salmon-1.0.0

echo "FINISHED at "$(date "+%Y-%m-%d-%H%M")""