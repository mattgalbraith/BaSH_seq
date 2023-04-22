#!/bin/bash -l
  
#SBATCH -o ./BWA_INDEX.out
#SBATCH -e ./BWA_INDEX.err
#SBATCH -D ./
#SBATCH -J BWA_INDEX
#SBATCH --partition=c2s30

SCRIPT_TITLE=BWA_INDEX.sh
SCRIPT_VERSION=0.1
# DATE: 11-21-2022
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script builds a bwa index from a genome fasta file (here = Genecode v33)
# Takes ~52 mins on Eureka c2s30 (up to 118 Gb RAM)

cd /home/matthew.galbraith/references/GRCh38_bwa

echo "Starting "$SCRIPT_TITLE" at "$(date "+%Y-%m-%d-%H%M")""

singularity run ~/Tools/Singularity/bwa-0.7.17-docker_amd64.sif bwa index \
	-p gencode.v33.basic \
	/home/matthew.galbraith/references/GRCh38/Gencode/genome/GRCh38.primary_assembly.genome.fa

echo "FINISHED at "$(date "+%Y-%m-%d-%H%M")""


# Usage:   bwa index [options] <in.fasta>

# Options: -a STR    BWT construction algorithm: bwtsw, is or rb2 [auto]
#          -p STR    prefix of the index [same as fasta name]
#          -b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]
#          -6        index files named as <in.fasta>.64.* instead of <in.fasta>.* 

# Warning: `-a bwtsw' does not work for short genomes, while `-a is' and
#          `-a div' do not work not for long genomes.