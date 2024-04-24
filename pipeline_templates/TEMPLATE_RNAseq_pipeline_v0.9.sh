#!/bin/bash

PIPELINE_TITLE=RNAseq_pipeline
PIPELINE_VERSION=0.9
# DATE: 06-30-2021
# AUTHOR: Matthew Galbraith
# SUMMARY: 
# This script is designed to be run once for each sample (incl. both reads for PE) and is executed by submitAll.sh once per sample with the appropriate arguments.
# Version 0.3+: Pipeline script is called from sbatch SUBMIT job; each stage is now also called by sbatch (with -W and --wrap options), rather than srun in earlier versions.
# Version 0.4+: Adds functions for catching SLURM stepd errors and checking SLURM exit codes
# Version 0.5: Activating stages 13 and 14 for HOMER makeTagDirectory and HOMER makeUCSCfile
# Version 0.6: Adding spike-in analysis
# Version 0.65: Adding SEQ_TYPE to stage 11 to trigger name sorting for PE bam
# Version 0.7: Use HISAT2 instead of TOPHAT2, step 7 is redunadant as bam is sorted in alignment step, but I didn't change much. Double sorting shouldn't hurt anything. Changed the output name in the alignment step (removed "sorted")
# Version 0.8_TEST: initial setup on HDC Eureka
# Version 0.8: initial working version for HDC Eureka 070621
# Version 0.9 102623: updating to run on Proton2 with Singularity/Apptainers
# 
# Steps to run pipeline:
# 1) In Project/ : create top-level working dirs: Project/raw_date and Project/analysis_date
# 2) In Project/analysis_date/ create sample_locations.txt (field 1 = SAMPLE_NAME; field 2 = path/to/raw_fastq_file.gz; field 3 = read1 / read2 labels)
# 3) Edit variables in top section of pipeline script and (if needed) save a project-specific version (specify path to this version as indicated below)
# 4) From Project/analysis_date/ run analysis_setup.sh as follows:
# sh path/to/analysis_setup.sh <FULL/PATH/TO/PIPELINE_SCRIPT> <START_AT_STAGE> <END_AT_STAGE>
# 5) Then from Project/analysis_date/scripts run submitAll_START_AT_STAGE-END_AT_STAGE.sh as follows (command can be copied from within submit script):
# sbatch submitAll_START_AT_STAGE-END_AT_STAGE.sh

##### ARGUMENTS PASSED FROM submitAll.sh - DO NOT EDIT THIS SECTION #####

EXPECTED_ARGS=6
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo "Usage: sh "$PIPELINE_TITLE"v"$PIPELINE_VERSION".sh <THIS_SAMPLE_NAME> <THIS_ANALYSIS_DIR> <RAW_DIR> <START_AT_STAGE> <END_AT_STAGE> <SUBMIT_LOG>
        This should normally be called by submitAll script
        ${red}ERROR - expecting "$EXPECTED_ARGS" ARGS but "$#" were provided:${NC}
        "$@"
        "
        exit 1
fi

# global pipeline variables from command line via submitAll.sh (do not edit)
THIS_SAMPLE_NAME=$1
THIS_ANALYSIS_DIR=$2
RAW_DIR=$3
START_AT_STAGE=$4
END_AT_STAGE=$5
SUBMIT_LOG=$6

##### VARIABLES IN THIS SECTION CAN BE EDITED #####
## global pipeline variables ##
        SEQ_TYPE=PE                         # <PE/SE> for paired-end or single end
        STRAND_TYPE=strand-specific-fwd     # <strand-specific-fwd/strand-specific-rev/unstranded>  # use fwd for Nugen/Tecan; use rev for Illumina TruSeq stranded
        READ_LENGTH=150                     # <50/75/150>  # used to trim n+1 read if present
        PIPELINE_TYPE=RNAseq
        SPIKE_IN=none                       # <none/both/dm/ercc> external RNA spike-in details; used to select correct index(s) and gtf(s) for alignment and counting - check that correct paths are ebtered for stages 4 and 11 below
        PIPELINE_SCRIPTS_DIR=/data1/matt_testing/BaSH_seq_test/BaSH_seq/RNAseq/stageScripts  # MAY NEED TO CUSTOMIZE
        QC_DIR_NAME="$THIS_ANALYSIS_DIR"/QC
        ALIGNMENT_DIRNAME="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Alignments   # NEED TO IMPLEMENT IN STAGE 4
        COUNTS_DIRNAME="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Counts
        TRACKS_DIRNAME="$THIS_ANALYSIS_DIR"/Tracks
        COMMENTS="PIPELINE RUN ON PROTON2"     # Add comments explaining this version of the analysis
        # NEW for 0.9: NEED TO GET RAW DIR FOR MOUNTING TO CONTAINERS THAT NEED TO READ SYMLINKS TO FASTQ FILES, eg FASTQC
        # RAW_DIR is passed to Pipeline script via Analysis_seetup.sh script
        # NEED TO ADD THIS TO STAGES THAT NEED TO SEE RAW_DIR
        # NEW for 0.9: Also need to mount references dir for stages that need to read references etc, eg FASTQ_SCREEN, HISAT2
        REFS_DIR=/data1/references/

## STAGE-SPECIFIC VARIABLES TO EDIT HERE ##
# 0 PRE_FASTQC
        # uncommon options are either defaults or set in 0_PRE_FASTQC.sh
        FASTQC_SIF=/data1/containers/fastqc0.11.9.sif

# 1 FASTQ_MCF
        # common options
        MIN_SEQ_LENGTH=30       # -l option; default = 19
        MIN_QUAL=10             # -q option; default = 10; should change for Novoseq data with binned q-scores
        MIN_ADAPTER=1.6           # -s option; Log2 scale for adapter minimum-length-match; default=2.2 ie ~4.5)
        MIN_PERC_OCCUR=0.25     # -t option; % occurance threshold before adapter clipping; default=0.25)
        MAX_PERC_DIFF=10        # -p option; Maximum adapter difference percentage; default=10)
        CONTAMINANTS_FASTA=/data1/references/Contaminants/contaminants.fa # need to bind this to container
        BBTOOLS_SIF=/data1/containers/bbtools39.01.sif
        EAUTILS_SIF=/data1/containers/eautils1.04.807.sif
        # uncommon options are either defaults or set in 1_FASTQ_MCF.sh
        # see also: https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqMcf.md

# 2 POST_FASTQC
        # uncommon options are either defaults or set in 2_POST_FASTQC.sh

# 3 FASTQ SCREEN
        # common options:
        # uncommon options are either defaults or set in 3_FASTQ_SCREEN.sh
        FASTQSCREEN_CONF=/data1/matt_testing/BaSH_seq_test/BaSH_seq/misc/fastq_screen.conf
        FASTQSCREEN_SIF=/data1/containers/fastqscreen0.15.2.sif

# # 4 MAPPING TOPHAT2_ALIGN
#         # common options
#         BOWTIE_INDEX=$SHARED/Kohl/references/hg19/gencode/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/
#         BOWTIE_INDEX_BOTH="$SHARED"/Matt/Refs/Bowtie2_UCSC_hg19_dm6_ERCC/hg19_dm6_ERCC.genome       # INDEX BASENAME FOR HUMAN + DROSOPHILA and ERCC SPIKE-IN
#         BOWTIE_INDEX_DM6=not_yet_implemented                                                        # INDEX BASENAME FOR HUMAN + DROSOPHILA
#         BOWTIE_INDEX_ERCC=not_yet_implemented                                                       # INDEX BASENAME FOR HUMAN + ERCC
#         # ^ use custom indexes when using Drosophila and/or ERCC spike-ins
#         # BASENAME FOR HUMAN: "$SHARED"/references/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
#         # BASENAME FOR MOUSE:
#         TOPHAT_GTF=$SHARED/Kohl/references/hg19/gencode/annotations/basic/gencode.v32lift37.basic.annotation.no_chr.gtf # renamed, from gencode
#         TOPHAT_GTF_BOTH="$SHARED"/Matt/Refs/Bowtie2_UCSC_hg19_dm6_ERCC/hg19_dm6_genes.gtf           # COMBINED GTF for hg19 and dm6; no need for ERCC as no introns
#         TOPHAT_GTF_DM6="$SHARED"/Matt/Refs/Bowtie2_UCSC_hg19_dm6_ERCC/hg19_dm6_genes.gtf            # same as for BOTH
#         # ^ use custom gtf when using Drosophila spike-ins
#         # FOR HUMAN: "$SHARED"/references/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf         # from igenomes
#         TRANSCRIPTOME_INDEX=$SHARED/Kohl/references/hg19/gencode/annotations/basic/transcriptome_data/known
#         TRANSCRIPTOME_INDEX_BOTH="$SHARED"/Matt/Refs/Bowtie2_UCSC_hg19_dm6_ERCC/transcriptome_index_bt2/hg19_dm6_known  # COMBINED hg19 + dm6
#         TRANSCRIPTOME_INDEX_DM6="$SHARED"/Matt/Refs/Bowtie2_UCSC_hg19_dm6_ERCC/transcriptome_index_bt2/hg19_dm6_known   # same as for BOTH
#         # ^ use custom index when using Drosophila spike-ins
#         # ^ created using this cmd: tophat -G /gpfs/jespinosa/shared/Matt/Refs/Bowtie2_UCSC_hg19_dm6_ERCC/hg19_dm6_genes.gtf --transcriptome-index=/gpfs/jespinosa/shared/Matt/Refs/Bowtie2_UCSC_hg19_dm6_ERCC/transcriptome_index_bt2/hg19_dm6_known /gpfs/jespinosa/shared/Matt/Refs/Bowtie2_UCSC_hg19_dm6_ERCC/hg19_dm6_ERCC.genome
#         # STD HUMAN: "$SHARED"/references/Homo_sapiens/UCSC/hg19/Annotation/Genes/transcriptome_index_bt2/genes # Should only be generated on first run
#         MATE_INNER_DIST=50      # Default is 50. Should be set based on library QC if possible. May require an initial alignment run to assess.
#         MATE_STD_DEV=20         # Default is 20.
#                                 # Leave set for SE data but will not be used.
#         # uncommon options are either defaults or set in 4_TOPHAT2_ALIGN.sh
        
# 4 MAPPING                             # Should test different insert limits and soft-clipping here - may also add second mapping stage...
        # common options
        ALIGNER="HISAT2_ALIGN"      # <BOWTIE2> or <HISTA2-DNA>
        ALIGNER_INDEX="$REFS_DIR"/GRCh38/Gencode/annotation/gencode.v33.basic # HUMAN
        ALIGNER_SS="$REFS_DIR"/GRCh38/Gencode/annotation/gencode.v33.basic.annotation.ss # HUMAN
        # ALIGNER_INDEX="$REFS_DIR"/GRCm38/Gencode/annotation/gencode.vM24.basic # MOUSE
        # ALIGNER_SS="$REFS_DIR"/GRCm38/Gencode/annotation/gencode.vM24.basic.annotation.ss # MOUSE
        # For Spike-ins:
        # ALIGNER_INDEX="$HOME"/references/Bowtie2_UCSC_hg19_ERCC_6/UCSC.hg19_ERCC_6.genome NOT CURRENTLY IMPLEMENTED
        # ^ this index includes the 6 ERCC spike-ins (3 labelled, 3 unlabelled)
        # For iGenomes ref (older analyses):
        # FOR HUMAN: "$REFS_DIR"/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
        # FOR MOUSE: 
        HISAT2_SIF=/data1/containers/hisat2_2.2.1.sif
        SAMTOOLS_SIF=/data1/containers/samtools1.16.1.sif
        # uncommon options are either defaults or set in 4_<ALIGNER>.sh

# 5 ADD RGID
       # common options
        PLATFORM=EXAMPLE_BaSH_seq
        DATE="$(date "+%m-%d-%Y")"
        PI=EXAMPLE_BaSH_seq
        LIBRARY="$SEQ_TYPE"_"$STRAND_TYPE"
        SEQ_CORE=EXAMPLE_BaSH_seq
        SEQ_ID=NA
        EXPERIMENT=stranded_RNAseq
        PICARD_SIF=/data1/containers/picard3.0.0.sif
        # uncommon options are either defaults or set in 5_
        PICARD_MEM_5=8G      # java memory option
        
# 6 FILTER
        MIN_MAPQ=10
        SAMTOOLS_SIF=/data1/containers/samtools1.16.1.sif
# 7 SORT
        PICARD_MEM_7=16G     # java memory option
                                                                        # need to check variables in script and move here
# 8 MARK OR REMOVE DUPLICATES
        PICARD_MEM_8=64G     # java memory options
        DUPLICATES=mark           # <mark> or <remove>
                                                                        # need to check variables in script and move here
# 9 ALIGNMENT METRICS
        PICARD_MEM_9=64G
        # REF_FILE="$HOME"/references/Picard/ucsc.hg19.fasta.gz # not currently used

# 10 RSEQC
        REFSEQ_BED=/data1/references/rseqc/hg38_Gencode_V33.bed # HUMAN
        # REFSEQ_BED=/data1/references/rseqc/mm10_RefSeq.bed # MOUSE
        HOUSEKEEPING_BED=/data1/references/rseqc/hg38.HouseKeepingGenes.bed # HUMAN
        # HOUSEKEEPING_BED=/data1/references/rseqc/mm10.HouseKeepingGenes.bed # MOUSE
        RSEQC_SIF=/data1/containers/rseqc5.0.1.sif
        RSEQC_REFS=/data1/references/rseqc


# 11 HTSeq COUNT
        GTF=/data1/references/GRCh38/Gencode/annotation/gencode.v33.basic.annotation.gtf # HUMAN
        # GTF=/data1/references/GRCm38/Gencode/annotation/gencode.vM24.basic.annotation.gtf # MOUSE
        #
        DM_GTF=/data1/references/dm6/Drosophila_melanogaster/UCSC/dm6/Annotation/Genes/genes.gtf           # from igenomes
        ERCC_GTF=/data1/references/ERCC/ERCC92/ERCC92.gtf                                                  # from https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip
        # ^ additional gtfs only used with Drosophila and/or ERCC spike-ins (do not comment out or args for stage get shifted)
        FEATURETYPE=exon
        IDATTR=gene_id # extracts ensemblID, which are unique
        HT_MINAQUAL=10
        MODE=intersection-nonempty
        ORDER=pos       # For paired-end data: <name> or <pos>; Default is name
        # BUFFER=       # --max-reads-in-buffer=<number>; when paired sorted by position; default is 30000000       # NOT IMPLEMENTED YET
        HTSEQ_SIF=/data1/containers/htseq2.0.2.sif

# NOT IN USE FOR RNAseq:
# # 12 HOMER MAKE TAG DIR
#         # Stages 11 and 12 should normally be run together to ensure same settings ie -fragLength, -tbp
#         HOMER_DIRNAME="$THIS_ANALYSIS_DIR"/HOMER
#         REF=hg19            # used for -genome <genome version> (required for -checkGC)
#         FRAG_LENGTH_s11=given  # "none" to use default; "given" to use read length (using custom BED format instead where read length = TLEN from BAM)
#                                # use given for PE data to ensure fragment lenghts are correct (= read length in custom BED file)
#                                ## HOMER does not fully support visualization for paired-end RNA-seq but if using fragLength given, it will plot to fist splice junction to preserve sharp exon boundaries. 
#         TAGS_PER_s11="all"      # used for -tbp INT  (default: no limit) number of tags per position to keep; "all" to use default
#                             # can also be set for makeUCSCfile
#         # -format {sam/bam/bed}     # can be used to force format if HOMER guesses incorrectly
#         # -d <tag directory> [tag directory 2] ... (add Tag directory to new tag directory)     # use to combine tag directories eg if plotting merged data
#         # -unique   (default) keep only uniquely aligned reads (secondary flag unset and MAPQ >10)
#         # -tbp INT  (default: no limit) number of tags per position to keep
#         # -single (Create a single tags.tsv file for all "chromosomes" - i.e. if >100 chromosomes)
#         # -checkGC  check Sequence bias, requires "-genome"
#         # If you use paired-end reads for ChIP-Seq or RNA-Seq, HOMER will treat each half of the read separately (and count each as 0.5 reads), which works well for a number of applications.
#         # If you are using stranded paired-end reads, make sure to specify "-sspe" so that HOMER will correctly interpret the intended strand for the 2nd read in the mate-pair.

# # 13 HOMER MAKE UCSC FILE
#         # Stages 11 and 12 should normally be run together to ensure same settings ie -fragLength, -tbp
#         FRAG_LENGTH_s12="$FRAG_LENGTH_s11"            # used for -fragLength <# | auto | given> option; default: auto; "none" to use default or should set to "$FRAG_LENGTH_s11" to match stage 11; may have to use given for PE RNAseq to prevent extension past splice sites (will give crisp borders at exons but not plot 3' parts of reads after junctions)
#         REF=hg19
#         TAGS_PER_s12="$TAGS_PER_s11"                 # used for -tbp <#> option; "all" to use default: no limit; if not filtering, should be set to "$TAGS_PER_s11" to match what was used for stage 11
#         HOMER_TAG_DIRNAME="$HOMER_DIRNAME"/"$THIS_SAMPLE_NAME"_TagDirectory_frag-"$FRAG_LENGTH_s11"_tbp-"$TAGS_PER_s11"
#         NORM_s12=1e6                    # used for -norm <#> option
#         RES_s12=1                       # used for -res <#>; resoultion in bp
#         ## -fsize <#>                               # hardcoded in script for now
#         # -normLength <#>                           # not implemented
#         ## -strand <both/separate>                  # using switch based on STRAND_TYPE
#         # Output format
#             # use switch for bedgraph (default for HOMER) vs bigwig     # not currently implemented


##### DO NOT EDIT PAST HERE FOR NORMAL USAGE #####

blue="\033[0;36m"
green="\033[0;32m"
red="\033[0;31m"
NC="\033[0m"    # no color (turn off color)

### Functions for error checking etc ###

function getJobID {
    # Gets job ID from submit script log file - but only after sbatch -W finishes
    # Probably could get from STAGE_OUTPUT file before job exits
    # Usage: getJobID <$JOBNAME>
    JOB_ID="$(grep -A 3 "Running $1" "$SUBMIT_LOG" | grep "Submitted batch job" | awk '{print $4}')" # This will need to be stage- and sample-specific
}

function slurmErrorCheck {
    # This looks at .err file after sbatch job has completed or failed to check for errors that do not pass non-zero exit codes to pipeline
    # Usage: slurmErrorCheck <$STAGE_ERROR> <$JOB_ID>
    if [ $(grep -c "slurmstepd: error:" "$1"."$2".*.err) -gt 0 ]
        then
            echo -e "${red}slurmstepd ERROR -check "$1"."$2".*.err ${NC}"
            grep "slurmstepd: error:" "$1"."$2".*.err
            OUTPUT_STATUS=1
    fi
}

function slurmExitCodeCheck {
    # This uses sacct to check exit codes for the job(s); there are 2 codes separated by ":" for the original job (called by called by sbatch --wrap inside pipeline script, the job.batch (called by sbatch --wrap inside pipeline script), and the job.0 (called by srun inside stage script)
    # Usage: slurmExitCodeCheck <$JOB_ID>
    if [ $(sacct -p -n -j $1 | awk -F'[|:]' '{print $7}' | grep -c -v "0") -gt 0 ] || [ $(sacct -p -n -j $1 | awk -F'[|:]' '{print $8}' | grep -c -v "0") -gt 0 ]
        then
            echo -e "${red}ERROR - non-zero exit code for job or job step:${NC}"
            sacct -j $1
            OUTPUT_STATUS=1
    fi
}

function copyOutErrLog {
    # Copy standard output and standard error to logfile
    # Usage: copyOutErrLog <Title> <$STAGE_ERROR> <$JOB_ID>
    echo -e "--------------------\nStage "$1"\n--------------------" >> "$LOGFILE"
    for i in $(ls "$2"."$3".*.out)
    do
        echo -ne "OUTPUT FILE: "$i"" >> "$LOGFILE"
        cat $i >> "$LOGFILE"
        rm $i
    done

    for i in $(ls "$2"."$3".*.err)
    do
        echo -e "\nERROR FILE: "$i"" >> "$LOGFILE"
        cat "$i" >> "$LOGFILE"
        rm "$i"
        echo -e "\n" >> "$LOGFILE"
    done
}

function checkOutputStatus {
    # Checks $OUTPUT_STATUS from sbatch command and exits if non-zero
    # Usage: checkOutputStatus <MESSAGE or $STAGE_NAME>
    if [ "$OUTPUT_STATUS" -ne 0 ]
        then
            echo -e "${red}Stage "$1" failed!!${NC}"
            OUTPUT_STATUS=""

            exit 1
    fi
}
###


echo -e "----------\n${blue}Starting "$PIPELINE_TITLE" v"$PIPELINE_VERSION" stage(s) "$START_AT_STAGE" to "$END_AT_STAGE" for "$THIS_SAMPLE_NAME"${NC}\n----------"



# Check for Sample FASTQ file(s) and set up FASTQR1 and/or FASTQR2 variables for PE/SE data
FASTQR1_FILE="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Raw/"$THIS_SAMPLE_NAME"_R1.fastq.gz
# for Eureka use basename of FASTQR1_FILE_ORIG instead of renaming; MAY NEED TO MOVE OR MAKE FUNCTION FOR ONLY STAGES THAT CARE ABOUT FASTQ
# FASTQR1_FILE=/tmp/$(basename "$THIS_ANALYSIS_DIR")/Sample_"$THIS_SAMPLE_NAME"/Raw/$(basename "$FASTQR1_FILE_ORIG")
if [ $START_AT_STAGE -lt 5 ] # Adding this to skip looking for FASTQ files after stage 4 (in case they have been moved to save disk space)
        then
        if [ ! -f "$FASTQR1_FILE" ]
        then
            echo -e "${red}FASTQR1 file does not exist: "$FASTQR1_FILE"
            "
            exit 1
        else
            echo "FASTQR1 file found: "$FASTQR1_FILE"
            "
        fi

        if [ $SEQ_TYPE = "SE" ]
        then
            echo "Running "$PIPELINE_TYPE" pipeline in single-end (SE) mode
            "
            FASTQR2_FILE="none"  # Set to "none" rather than empty to avoid sending empty variable (eg to FASTQ_MCF / FASTQ_SCREEN / BOWTIE2) # Should not matter if using $SEQ_TYPE switch

        elif [ $SEQ_TYPE = "PE" ]
        then
            echo "Running "$PIPELINE_TYPE" pipeline in paired-end (PE) mode
            "
            FASTQR2_FILE="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Raw/"$THIS_SAMPLE_NAME"_R2.fastq.gz

            if [ ! -f "$FASTQR2_FILE" ]
            then
                echo -e "${red}FASTQR2 file does not exist: "$FASTQR2_FILE"
                "
                exit 1
            else
                echo "FASTQR2 file found: "$FASTQR2_FILE"
                "
            fi
        fi
fi

##### SETTING UP DIRECTORIES #####
# Create std_out_err directory
if [ ! -d "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/ ]
then
        echo "making std_out_err directory "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/
        "
        mkdir "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/
fi


##### Capture snapshot of stage scripts to be used for this pipeline run #####
# This will ensure that old versions are preserved in case of major changes
SCRIPTS_BAK_DIR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/stageScripts_bak/"$PIPELINE_TYPE"-"$(basename "$THIS_ANALYSIS_DIR")"-Stage"$START_AT_STAGE"-"$END_AT_STAGE"-"$THIS_SAMPLE_NAME"-"$(date "+%Y-%m-%d-%H%M")"
mkdir --parents "$SCRIPTS_BAK_DIR"
for i in $(seq ${START_AT_STAGE} 1 ${END_AT_STAGE})
do
    cp "$PIPELINE_SCRIPTS_DIR"/${i}_*.sh  "$SCRIPTS_BAK_DIR"
done

##### PIPELINE LOGFILE #####
LOGFILE="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$PIPELINE_TYPE"-"$(basename "$THIS_ANALYSIS_DIR")"-Stage"$START_AT_STAGE"-"$END_AT_STAGE"-"$THIS_SAMPLE_NAME"-"$(date "+%Y-%m-%d-%H%M")".log
echo "
#######################################
"$PIPELINE_TITLE" v"$PIPELINE_VERSION"
Stage(s): "$START_AT_STAGE"-"$END_AT_STAGE"
Run: $(date)
Sample: "$THIS_SAMPLE_NAME"
#######################################
COMMENTS:
"$COMMENTS"

GLOBAL SETTINGS:
From command line (via analysis_setup.sh):
Sample name: "$THIS_SAMPLE_NAME"
Analysis directory: "$THIS_ANALYSIS_DIR"
Raw file directory: "$RAW_DIR"
User account: "$THIS_USER_ACCOUNT"
Begin at stage: "$START_AT_STAGE"
Finish at stage: "$END_AT_STAGE"
Set within "$PIPELINE_TITLE"v"$PIPELINE_VERSION".sh
Sequencing type: "$SEQ_TYPE"
Strand type: "$STRAND_TYPE"
Read length: "$READ_LENGTH"
Pipeline scripts directory: "$PIPELINE_SCRIPTS_DIR"
QC directory: "$QC_DIR_NAME"
Tracks directory: "$TRACKS_DIRNAME"
Read 1 fastq file: "$FASTQR1_FILE"
Read 2 fastq file: "$FASTQR2_FILE"
" > "$LOGFILE"
# ADD some global report params? Comments?

##### PIPELINE STAGES #####

# 0) PRE_FASTQC
        # run the stage
        STAGE_NAME="0-PRE_FASTQ"
        
        if [ $START_AT_STAGE -eq 0 ]
        then
                echo -e "${blue}Starting stage "$STAGE_NAME" for "$THIS_SAMPLE_NAME"${NC}
                "

                if [ ! -d "$QC_DIR_NAME" ]
                then
                      echo "making QC directory: "$QC_DIR_NAME"
                      "
                      mkdir --parents "$QC_DIR_NAME"
                fi

                #sh $PIPELINE_SCRIPTS_DIR/0_PRE_FASTQC.sh
                #Usage: 0_PRE_FASTQC.sh <SAMPLE_DIR> <SAMPLE_NAME> <FASTQR1_FILE/FASTQR2_FILE> <QC_DIR_NAME> <OUT_DIR_NAME> <THREADS> <FASTQC_SIF> <THIS_ANALYSIS_DIR> <RAW_DIR> 

                # Read 1
                JOB_NAME=stage-"$STAGE_NAME"_R1-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                echo -e "${blue}Running "$JOB_NAME"${NC}
                "
                sbatch -W \
                        --job-name="$JOB_NAME" \
                        --output="$STAGE_OUTPUT".%j.%N.out \
                        --error="$STAGE_ERROR".%j.%N.err \
                        --partition=defq \
                        --time=10:00:00 \
                        --nodes=1 \
                        --ntasks=1 \
                        --cpus-per-task=8 \
                        --mem-per-cpu=4G \
                        --wrap="\
                                sh "$PIPELINE_SCRIPTS_DIR"/0_PRE_FASTQC.sh "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ "$THIS_SAMPLE_NAME" "$FASTQR1_FILE" "$QC_DIR_NAME" FASTQC_pre_filtered 8 "$FASTQC_SIF" "$THIS_ANALYSIS_DIR" "$RAW_DIR"\
                                "

                # Catch output status
                OUTPUT_STATUS=$?

                # Get Job ID
                getJobID "$JOB_NAME"

                # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                # Check for exit code error: <function> <$JOB_ID>
                slurmExitCodeCheck "$JOB_ID"
                
                # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                copyOutErrLog "$STAGE_NAME"_read1 "$STAGE_ERROR" "$JOB_ID"

                # check output status: <function> <stage label>
                checkOutputStatus "$STAGE_NAME"_read1
                

                # Check for PE and run Read 2
                if [ $SEQ_TYPE = "PE" ]
                then
                    JOB_NAME=stage-"$STAGE_NAME"_R2-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                    STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                    STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                    echo -e "${blue}Running "$JOB_NAME"${NC}
                    "
                    sbatch -W \
                            --job-name="$JOB_NAME" \
                            --output="$STAGE_OUTPUT".%j.%N.out \
                            --error="$STAGE_ERROR".%j.%N.err \
                            --partition=defq \
                            --time=10:00:00 \
                            --nodes=1 \
                            --ntasks=1 \
                            --cpus-per-task=8 \
                            --mem-per-cpu=4G \
                            --wrap="\
                                    sh "$PIPELINE_SCRIPTS_DIR"/0_PRE_FASTQC.sh "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ "$THIS_SAMPLE_NAME" "$FASTQR2_FILE" "$QC_DIR_NAME" FASTQC_pre_filtered 8 "$FASTQC_SIF" "$THIS_ANALYSIS_DIR" "$RAW_DIR"\
                                    "

                    # Catch output status
                    OUTPUT_STATUS=$?

                    # Get Job ID
                    getJobID "$JOB_NAME"

                    # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                    slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                    # Check for exit code error: <function> <$JOB_ID>
                    slurmExitCodeCheck "$JOB_ID"
                    
                    # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                    copyOutErrLog "$STAGE_NAME"_read2 "$STAGE_ERROR" "$JOB_ID"

                    # check output status: <function> <stage label>
                    checkOutputStatus "$STAGE_NAME"_read2
                
                fi
                
                echo -e "${green}Stage "$STAGE_NAME" complete for "$THIS_SAMPLE_NAME"${NC}
                "
                # update the stage
                START_AT_STAGE=$(($START_AT_STAGE+1))
                STAGE_NAME=""
                STAGE_OUTPUT=""
                STAGE_ERROR=""
        fi

        # check if we need to exit
        if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
        then
                echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
                exit 0
        fi


# 1) FASTQ-MCF
        # run the stage
        STAGE_NAME="1-FASTQ-MCF"
        if [ $START_AT_STAGE -eq 1 ]
        then
                JOB_NAME=stage-"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                echo -e "${blue}Running "$JOB_NAME"${NC}
                "

                # make a Processed directory inside the sample
                if [ ! -d "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Processed/ ]
                then
                        mkdir "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Processed/
                fi

                #sh $PIPELINE_SCRIPTS_DIR/1_FASTQ_MCF.sh
                #Usage: 1_FASTQ_MCF.sh <SEQ_TYPE> <READ_LENGTH> <SAMPLE_DIR> <SAMPLE_NAME> <FASTQR1_FILE> <FASTQR2_FILE> <QC_DIR_NAME> <MIN_QUAL> <MIN_SEQ_LENGTH> <MIN_ADAPTER> <MIN_PERC_OCCUR> <MAX_PERC_DIFF> <CONTAMINANTS_FASTA> <BBTOOLS_SIF> <EAUTILS_SIF> <THIS_ANALYSIS_DIR> <RAW_DIR>                            

                sbatch -W \
                        --account="$THIS_USER_ACCOUNT" \
                        --job-name="$JOB_NAME" \
                        --output="$STAGE_OUTPUT".%j.%N.out \
                        --error="$STAGE_ERROR".%j.%N.err \
                        --partition=defq \
                        --time=10:00:00 \
                        --nodes=1 \
                        --ntasks=1 \
                        --cpus-per-task=8 \
                        --mem-per-cpu=10G \
                        --wrap="\
                                sh "$PIPELINE_SCRIPTS_DIR"/1_FASTQ_MCF_BB.sh \
                                "$SEQ_TYPE" \
                                "$READ_LENGTH" \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ \
                                "$THIS_SAMPLE_NAME" \
                                "$FASTQR1_FILE" \
                                "$FASTQR2_FILE" \
                                "$QC_DIR_NAME" \
                                "$MIN_QUAL" \
                                "$MIN_SEQ_LENGTH" \
                                "$MIN_ADAPTER" \
                                "$MIN_PERC_OCCUR" \
                                "$MAX_PERC_DIFF" \
                                "$CONTAMINANTS_FASTA" \
                                "$BBTOOLS_SIF" \
                                "$EAUTILS_SIF" \
                                "$THIS_ANALYSIS_DIR" \
                                "$RAW_DIR"
                                "

                # Catch output status
                OUTPUT_STATUS=$?

                # Get Job ID
                getJobID "$JOB_NAME"

                # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                # Check for exit code error: <function> <$JOB_ID>
                slurmExitCodeCheck "$JOB_ID"
                
                # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"

                # check output status: <function> <stage label>
                checkOutputStatus "$STAGE_NAME"
       
                echo -e "${green}"$JOB_NAME" complete.${NC}
                "
                # update the stage
                START_AT_STAGE=$(($START_AT_STAGE+1))
                STAGE_NAME=""
                STAGE_OUTPUT=""
                STAGE_ERROR=""
        fi
        
        # check if we need to exit
        if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
        then
                echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
                exit 0
        fi


# 2) POST_FASTQC 
        # run the stage
        STAGE_NAME="2-POST_FASTQC"
        if [ $START_AT_STAGE -eq 2 ]
        then
                echo -e "${blue}Starting stage "$STAGE_NAME" for "$THIS_SAMPLE_NAME"${NC}
                "

                #sh $PIPELINE_SCRIPTS_DIR/2_POST_FASTQC.sh
                #Usage: 2_POST_FASTQC.sh <SAMPLE_DIR> <SAMPLE_NAME> <FASTQ_FILE_NAME> <QC_DIR_NAME> <OUT_DIR_NAME> <THREADS> <FASTQC_SIF> <THIS_ANALYSIS_DIR> <RAW_DIR>

                # Read 1
                JOB_NAME=stage-"$STAGE_NAME"_R1-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                echo -e "${blue}Running "$JOB_NAME"${NC}
                "
                sbatch -W \
                        --account="$THIS_USER_ACCOUNT" \
                        --job-name="$JOB_NAME" \
                        --output="$STAGE_OUTPUT".%j.%N.out \
                        --error="$STAGE_ERROR".%j.%N.err \
                        --partition=defq \
                        --time=10:00:00 \
                        --nodes=1 \
                        --ntasks=1 \
                        --cpus-per-task=8 \
                        --mem-per-cpu=4G \
                        --wrap="\
                                sh "$PIPELINE_SCRIPTS_DIR"/2_POST_FASTQC.sh "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ "$THIS_SAMPLE_NAME" "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Processed/trimmed_"$(basename "$FASTQR1_FILE")" "$QC_DIR_NAME" FASTQC_post_filtered 8 "$FASTQC_SIF" "$THIS_ANALYSIS_DIR" "$RAW_DIR"\
                                "

                # Catch output status
                OUTPUT_STATUS=$?

                # Get Job ID
                getJobID "$JOB_NAME"

                # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                # Check for exit code error: <function> <$JOB_ID>
                slurmExitCodeCheck "$JOB_ID"
                
                # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                copyOutErrLog "$STAGE_NAME"_read1 "$STAGE_ERROR" "$JOB_ID"

                # check output status: <function> <stage label>
                checkOutputStatus "$STAGE_NAME"_read1


                # Check for PE and run Read 2
                if [ $SEQ_TYPE = "PE" ]
                then
                    JOB_NAME=stage-"$STAGE_NAME"_R2-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                    STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                    STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                    echo -e "${blue}Running "$JOB_NAME"${NC}
                    "
                    sbatch -W \
                            --account="$THIS_USER_ACCOUNT" \
                            --job-name="$JOB_NAME" \
                            --output="$STAGE_OUTPUT".%j.%N.out \
                            --error="$STAGE_ERROR".%j.%N.err \
                            --partition=defq \
                            --time=10:00:00 \
                            --nodes=1 \
                            --ntasks=1 \
                            --cpus-per-task=8 \
                            --mem-per-cpu=4G \
                            --wrap="\
                                    sh "$PIPELINE_SCRIPTS_DIR"/2_POST_FASTQC.sh "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ "$THIS_SAMPLE_NAME" "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Processed/trimmed_"$(basename "$FASTQR2_FILE")" "$QC_DIR_NAME" FASTQC_post_filtered 8 "$FASTQC_SIF" "$THIS_ANALYSIS_DIR" "$RAW_DIR"\
                                    "

                    # Catch output status
                    OUTPUT_STATUS=$?

                    # Get Job ID
                    getJobID "$JOB_NAME"

                    # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                    slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                    # Check for exit code error: <function> <$JOB_ID>
                    slurmExitCodeCheck "$JOB_ID"
                    
                    # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                    copyOutErrLog "$STAGE_NAME"_read2 "$STAGE_ERROR" "$JOB_ID"

                    # check output status: <function> <stage label>
                    checkOutputStatus "$STAGE_NAME"_read2

                fi

                echo -e "${green}Stage "$STAGE_NAME" complete for "$THIS_SAMPLE_NAME"${NC}
                "
                # update the stage
                START_AT_STAGE=$(($START_AT_STAGE+1))
                STAGE_NAME=""
                STAGE_OUTPUT=""
                STAGE_ERROR=""
        fi
        
        # check if we need to exit
        if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
        then
                echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
                exit 0
        fi


# 3) FASTQ SCREEN
        # run the stage
        STAGE_NAME="3-FASTQ_SCREEN"
        if [ $START_AT_STAGE -eq 3 ]
        then
                JOB_NAME=stage-"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                echo -e "${blue}Running "$JOB_NAME"${NC}
                "

                #sh $PIPELINE_SCRIPTS_DIR/3_FASTQ_SCREEN.sh
                #Usage: 3_FASTQ_SCREEN.sh <SEQ_TYPE> <SAMPLE_DIR> <SAMPLE_NAME> <FASTQR1_FILE> <FASTQR2_FILE> <QC_DIR_NAME> <THREADS> <FASTQSCREEN_CONF> <FASTQSCREEN_SIF> <THIS_ANALYSIS_DIR> <RAW_DIR> <REFS_DIR>
                               
                sbatch -W \
                        --account="$THIS_USER_ACCOUNT" \
                        --job-name="$JOB_NAME" \
                        --output="$STAGE_OUTPUT".%j.%N.out \
                        --error="$STAGE_ERROR".%j.%N.err \
                        --partition=defq \
                        --time=10:00:00 \
                        --nodes=1 \
                        --ntasks=1 \
                        --cpus-per-task=8 \
                        --mem-per-cpu=4G \
                        --wrap="\
                                sh "$PIPELINE_SCRIPTS_DIR"/3_FASTQ_SCREEN.sh \
                                "$SEQ_TYPE" \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ \
                                "$THIS_SAMPLE_NAME" \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Processed/trimmed_"$(basename "$FASTQR1_FILE")" \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Processed/trimmed_"$(basename "$FASTQR2_FILE")" \
                                "$QC_DIR_NAME" \
                                16 \
                                "$FASTQSCREEN_CONF" \
                                "$FASTQSCREEN_SIF" \
                                "$THIS_ANALYSIS_DIR" \
                                "$RAW_DIR" \
                                "$REFS_DIR"\
                                "

                # Catch output status
                OUTPUT_STATUS=$?

                # Get Job ID
                getJobID "$JOB_NAME"

                # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                # Check for exit code error: <function> <$JOB_ID>
                slurmExitCodeCheck "$JOB_ID"
                
                # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"

                # check output status: <function> <stage label>
                checkOutputStatus "$STAGE_NAME"
       
                echo -e "${green}"$JOB_NAME" complete.${NC}
                "
                # update the stage
                START_AT_STAGE=$(($START_AT_STAGE+1))
                STAGE_NAME=""
                STAGE_OUTPUT=""
                STAGE_ERROR=""
        fi
        
        # check if we need to exit
        if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
        then
                echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
                exit 0
        fi


# # 4) MAPPING using TOPHAT2
#         # run the stage
#         STAGE_NAME="4-TOPHAT2-ALIGN"
#         if [ $START_AT_STAGE -eq 4 ]
#         then
#                 JOB_NAME=stage-"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
#                 STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
#                 STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
#                 echo -e "${blue}Running "$JOB_NAME"${NC}
#                 "
# 
#                 # switch to select correct index and GTF
#                 if [ $SPIKE_IN = none ]
#                 then
#                     INDEX="$BOWTIE_INDEX"
#                     GTF="$TOPHAT_GTF"
#                     TX_INDEX="$TRANSCRIPTOME_INDEX"
#                 elif [ $SPIKE_IN = both ]
#                 then
#                     INDEX="$BOWTIE_INDEX_BOTH"
#                     GTF="$TOPHAT_GTF_BOTH"
#                     TX_INDEX="$TRANSCRIPTOME_INDEX_BOTH"
#                 elif [ $SPIKE_IN = dm ]
#                 then
#                     INDEX="$BOWTIE_INDEX_DM6"
#                     GTF="$TOPHAT_GTF_DM6"
#                     TX_INDEX="$TRANSCRIPTOME_INDEX_DM6"
#                 elif [ $SPIKE_IN = ercc ]
#                 then
#                     INDEX="$BOWTIE_INDEX_ERCC"
#                     GTF="$TOPHAT_GTF"
#                     TX_INDEX="$TRANSCRIPTOME_INDEX"
#                 fi
# 
#                 #sh $PIPELINE_SCRIPTS_DIR/4_TOPHAT2_ALIGN.sh
#                 #Usage: 4_TOPHAT2_ALIGN.sh <SEQ_TYPE> <SAMPLE_DIR> <SAMPLE_NAME> <FASTQR1_FILE> <FASTQR2_FILE> <INDEX> <GTF> <TX_INDEX> <STRAND_TYPE> <MATE_INNER_DIST> <MATE_STD_DEV> <BAM_OUT_FILENAME> <OUT_DIR_NAME> <THREADS>
#                                            
#                 sbatch -W \
#                         --account="$THIS_USER_ACCOUNT" \
#                         --job-name="$JOB_NAME" \
#                         --output="$STAGE_OUTPUT".%j.%N.out \
#                         --error="$STAGE_ERROR".%j.%N.err \
#                         --partition=bigmem \
#                         --time=10:00:00 \
#                         --mem=400G \
#                         --nodes=1 \
#                         --cpus-per-task=12 \
#                         --ntasks=1 \
#                         --wrap="\
#                                 sh $PIPELINE_SCRIPTS_DIR/4_TOPHAT2_ALIGN.sh "$SEQ_TYPE" "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ "$THIS_SAMPLE_NAME" "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Processed/trimmed_"$(basename "$FASTQR1_FILE")" "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Processed/trimmed_"$(basename "$FASTQR2_FILE")" "$INDEX" "$GTF" "$TX_INDEX" "$STRAND_TYPE" "$MATE_INNER_DIST" "$MATE_STD_DEV" "$THIS_SAMPLE_NAME".mapped.no-rgid.bam Alignments 12\
#                                 "
# 
#                 # Catch output status
#                 OUTPUT_STATUS=$?
# 
#                 # Get Job ID
#                 getJobID "$JOB_NAME"
# 
#                 # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
#                 slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"
# 
#                 # Check for exit code error: <function> <$JOB_ID>
#                 slurmExitCodeCheck "$JOB_ID"
#                 
#                 # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
#                 copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"
# 
#                 # check output status: <function> <stage label>
#                 checkOutputStatus "$STAGE_NAME"
#        
#                 echo -e "${green}"$JOB_NAME" complete.${NC}
#                 "
#                 # update the stage
#                 START_AT_STAGE=$(($START_AT_STAGE+1))
#                 STAGE_NAME=""
#                 STAGE_OUTPUT=""
#                 STAGE_ERROR=""
#         fi
#         
#         # check if we need to exit
#         if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
#         then
#                 echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
#                 exit 0
#         fi


# 4) MAPPING using HISAT2
        # run the stage
        STAGE_NAME="4-"$ALIGNER""
        if [ $START_AT_STAGE -eq 4 ]
        then
                JOB_NAME=stage-"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                echo -e "${blue}Running "$JOB_NAME"${NC}
                "

                #sh $PIPELINE_SCRIPTS_DIR/4_HISAT2_ALIGN.sh
                #Usage: 4_HISAT2.sh <SEQ_TYPE> <STRAND_TYPE> <SAMPLE_DIR> <SAMPLE_NAME> <FASTQR1_FILE> <FASTQR2_FILE> <ALIGNER_INDEX> <ALIGNER_SS> <BAM_OUT_FILENAME> <OUT_DIR_NAME> <THREADS> <HISAT2_SIF> <SAMTOOLS_SIF> <THIS_ANALYSIS_DIR> <RAW_DIR> <REFS_DIR>
                                           
                sbatch -W \
                        --account="$THIS_USER_ACCOUNT" \
                        --job-name="$JOB_NAME" \
                        --output="$STAGE_OUTPUT".%j.%N.out \
                        --error="$STAGE_ERROR".%j.%N.err \
                        --partition=defq \
                        --time=10:00:00 \
                        --nodes=1 \
                        --ntasks=1 \
                        --cpus-per-task=10 \
                        --mem-per-cpu=32G \
                        --wrap="\
                                sh $PIPELINE_SCRIPTS_DIR/4_"$ALIGNER".sh \
                                "$SEQ_TYPE" \
                                "$STRAND_TYPE" \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ \
                                "$THIS_SAMPLE_NAME" \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Processed/trimmed_"$(basename "$FASTQR1_FILE")" \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Processed/trimmed_"$(basename "$FASTQR2_FILE")" \
                                "$ALIGNER_INDEX" \
                                "$ALIGNER_SS" \
                                "$THIS_SAMPLE_NAME".mapped.no-rgid.bam \
                                Alignments \
                                10 \
                                "$HISAT2_SIF" \
                                "$SAMTOOLS_SIF" \
                                "$THIS_ANALYSIS_DIR" \
                                "$RAW_DIR" \
                                "$REFS_DIR"\
                                "

                # Catch output status
                OUTPUT_STATUS=$?

                # Get Job ID
                getJobID "$JOB_NAME"

                # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                # Check for exit code error: <function> <$JOB_ID>
                slurmExitCodeCheck "$JOB_ID"

                # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"

                # check output status: <function> <stage label>
                checkOutputStatus "$STAGE_NAME"

                echo -e "${green}"$JOB_NAME" complete.${NC}
                "
                # update the stage
                START_AT_STAGE=$(($START_AT_STAGE+1))
                STAGE_NAME=""
                STAGE_OUTPUT=""
                STAGE_ERROR=""
        fi

        # check if we need to exit
        if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
        then
                echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
                exit 0
        fi
        


# 5) ADD RGID
        # run the stage
        STAGE_NAME="5-ADD_RGID"
        if [ $START_AT_STAGE -eq 5 ]
        then
                JOB_NAME=stage-"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                echo -e "${blue}Running "$JOB_NAME"${NC}
                "

                #sh $PIPELINE_SCRIPTS_DIR/5_ADD_RGID.sh
                #Usage: 5_ADD_RGID.sh <SAMPLE_NAME> <PLATFORM> <DATE> <PI (investigator code)> <LIBRARY (PE/SE)> <SEQ_CORE> <SEQ_ID> <EXPERIMENT> <SAMPLE_DIR> <ALIGNMENT_DIRNAME> <BAM_IN_FILENAME> <BAM_OUT_FILENAME> <PICARD_MEM> <PICARD_SIF> <THIS_ANALYSIS_DIR>
                              
                sbatch -W \
                        --account="$THIS_USER_ACCOUNT" \
                        --job-name="$JOB_NAME" \
                        --output="$STAGE_OUTPUT".%j.%N.out \
                        --error="$STAGE_ERROR".%j.%N.err \
                        --partition=defq \
                        --time=10:00:00 \
                        --nodes=1 \
                        --ntasks=1 \
                        --cpus-per-task=1 \
                        --mem-per-cpu=8G \
                        --wrap="\
                                sh $PIPELINE_SCRIPTS_DIR/5_ADD_RGID.sh \
                                "$THIS_SAMPLE_NAME" \
                                "$PLATFORM" "$DATE" "$PI" "$LIBRARY" "$SEQ_CORE" "$SEQ_ID" "$EXPERIMENT" \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ Alignments "$THIS_SAMPLE_NAME".mapped.no-rgid.bam "$THIS_SAMPLE_NAME".mapped.rgid.bam \
                                "$PICARD_MEM_5" \
                                "$PICARD_SIF" \
                                "$THIS_ANALYSIS_DIR"\
                                "

                # Catch output status
                OUTPUT_STATUS=$?

                # Get Job ID
                getJobID "$JOB_NAME"

                # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                # Check for exit code error: <function> <$JOB_ID>
                slurmExitCodeCheck "$JOB_ID"
                
                # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"

                # check output status: <function> <stage label>
                checkOutputStatus "$STAGE_NAME"
       
                echo -e "${green}"$JOB_NAME" complete.${NC}
                "
                # update the stage
                START_AT_STAGE=$(($START_AT_STAGE+1))
                STAGE_NAME=""
                STAGE_OUTPUT=""
                STAGE_ERROR=""
        fi

        # check if we need to exit
        if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
        then
                echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
                exit 0
        fi


# 6) MAPQ FILTER
        # run the stage
        STAGE_NAME="6-MAPQ_FILTER"
        if [ $START_AT_STAGE -eq 6 ]
        then
                JOB_NAME=stage-"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                echo -e "${blue}Running "$JOB_NAME"${NC}
                "

                #sh $PIPELINE_SCRIPTS_DIR/6_MAPQ_FILTER.sh
                #Usage: 6_MAPQ_FILTER.sh <SAMPLE_DIR> <ALIGNMENT_DIRNAME> <BAM_IN_FILENAME> <BAM_OUT_FILENAME> <MIN_MAPQ> <SAMTOOLS_SIF> <THIS_ANALYSIS_DIR>
                                
                sbatch -W \
                        --account="$THIS_USER_ACCOUNT" \
                        --job-name="$JOB_NAME" \
                        --output="$STAGE_OUTPUT".%j.%N.out \
                        --error="$STAGE_ERROR".%j.%N.err \
                        --partition=defq \
                        --time=10:00:00 \
                        --nodes=1 \
                        --ntasks=1 \
                        --cpus-per-task=1 \
                        --mem-per-cpu=8G \
                        --wrap="\
                                sh $PIPELINE_SCRIPTS_DIR/6_MAPQ_FILTER.sh \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ Alignments "$THIS_SAMPLE_NAME".mapped.rgid.bam "$THIS_SAMPLE_NAME".mapped.rgid.filtered.bam \
                                "$MIN_MAPQ" \
                                "$SAMTOOLS_SIF" \
                                "$THIS_ANALYSIS_DIR"\
                                "

                # Catch output status
                OUTPUT_STATUS=$?

                # Get Job ID
                getJobID "$JOB_NAME"

                # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                # Check for exit code error: <function> <$JOB_ID>
                slurmExitCodeCheck "$JOB_ID"
                
                # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"

                # check output status: <function> <stage label>
                checkOutputStatus "$STAGE_NAME"
       
                echo -e "${green}"$JOB_NAME" complete.${NC}
                "
                # update the stage
                START_AT_STAGE=$(($START_AT_STAGE+1))
                STAGE_NAME=""
                STAGE_OUTPUT=""
                STAGE_ERROR=""
        fi

        # check if we need to exit
        if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
        then
                echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
                exit 0
        fi


# 7) SORT BAM
        # run the stage
        STAGE_NAME="7-SORT_BAM"
        if [ $START_AT_STAGE -eq 7 ]
        then
                JOB_NAME=stage-"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                echo -e "${blue}Running "$JOB_NAME"${NC}
                "

                #sh $PIPELINE_SCRIPTS_DIR/7_SORT_BAM.sh
                #Usage: 7_SORT_BAM.sh <SAMPLE_DIR> <ALIGNMENT_DIRNAME> <BAM_IN_FILENAME> <BAM_OUT_FILENAME> <PICARD_MEM> <PICARD_SIF> <THIS_ANALYSIS_DIR>
                
                sbatch -W \
                        --account="$THIS_USER_ACCOUNT" \
                        --job-name="$JOB_NAME" \
                        --output="$STAGE_OUTPUT".%j.%N.out \
                        --error="$STAGE_ERROR".%j.%N.err \
                        --partition=defq \
                        --time=10:00:00 \
                        --nodes=1 \
                        --ntasks=1 \
                        --cpus-per-task=8 \
                        --mem-per-cpu=4G \
                        --wrap="\
                                sh "$PIPELINE_SCRIPTS_DIR"/7_SORT_BAM.sh \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ Alignments "$THIS_SAMPLE_NAME".mapped.rgid.filtered.bam "$THIS_SAMPLE_NAME".mapped.rgid.filtered.sorted.bam \
                                "$PICARD_MEM_7" \
                                "$PICARD_SIF" \
                                "$THIS_ANALYSIS_DIR"\
                                "

                # Catch output status
                OUTPUT_STATUS=$?

                # Get Job ID
                getJobID "$JOB_NAME"

                # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                # Check for exit code error: <function> <$JOB_ID>
                slurmExitCodeCheck "$JOB_ID"
                
                # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"

                # check output status: <function> <stage label>
                checkOutputStatus "$STAGE_NAME"
       
                echo -e "${green}"$JOB_NAME" complete.${NC}
                "
                # update the stage
                START_AT_STAGE=$(($START_AT_STAGE+1))
                STAGE_NAME=""
                STAGE_OUTPUT=""
                STAGE_ERROR=""
        fi

        # check if we need to exit
        if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
        then
                echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
                exit 0
        fi


# 8) MARK DUPLICATES
        # run the stage
        STAGE_NAME="8-MARK_DUPLICATES"
        if [ $START_AT_STAGE -eq 8 ]
        then
                JOB_NAME=stage-"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                echo -e "${blue}Running "$JOB_NAME"${NC}
                "

                #sh $PIPELINE_SCRIPTS_DIR/8_MARK_DUPLICATES.sh 
                #Usage: 8_MARK_DUPLICATES.sh <SAMPLE_DIR> <ALIGNMENT_DIRNAME> <BAM_IN_FILENAME> <BAM_OUT_FILENAME> <PICARD_MEM> <mark/remove duplicates> <PICARD_SIF> <THIS_ANALYSIS_DIR>

                sbatch -W \
                        --account="$THIS_USER_ACCOUNT" \
                        --job-name="$JOB_NAME" \
                        --output="$STAGE_OUTPUT".%j.%N.out \
                        --error="$STAGE_ERROR".%j.%N.err \
                        --partition=defq \
                        --time=10:00:00 \
                        --nodes=1 \
                        --ntasks=1 \
                        --cpus-per-task=8 \
                        --mem-per-cpu=16G \
                        --wrap="\
                                sh "$PIPELINE_SCRIPTS_DIR"/8_MARK_DUPLICATES.sh \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ Alignments "$THIS_SAMPLE_NAME".mapped.rgid.filtered.sorted.bam "$THIS_SAMPLE_NAME".mapped.rgid.filtered.sorted.dups_mark.bam \
                                "$PICARD_MEM_8" \
                                "$DUPLICATES" \
                                "$PICARD_SIF" \
                                "$THIS_ANALYSIS_DIR"\
                                "

                # Catch output status
                OUTPUT_STATUS=$?

                # Get Job ID
                getJobID "$JOB_NAME"

                # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                # Check for exit code error: <function> <$JOB_ID>
                slurmExitCodeCheck "$JOB_ID"
                
                # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"

                # check output status: <function> <stage label>
                checkOutputStatus "$STAGE_NAME"
       
                echo -e "${green}"$JOB_NAME" complete.${NC}
                "
                # update the stage
                START_AT_STAGE=$(($START_AT_STAGE+1))
                STAGE_NAME=""
                STAGE_OUTPUT=""
                STAGE_ERROR=""
        fi

        # check if we need to exit
        if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
        then
                echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
                exit 0
        fi


# 9) ALIGNMENT METRICS
        # run the stage
        STAGE_NAME="9-ALIGNMENT_METRICS"
        if [ $START_AT_STAGE -eq 9 ]

        then
                JOB_NAME=stage-"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                echo -e "${blue}Running "$JOB_NAME"${NC}
                "

                if [ "$SEQ_TYPE" = "PE" ]
                    then
                        echo "Collecting paired-end alignment metrics"
                        ALIGNMENT_METRICS_SCRIPT="9_PE_ALIGNMENT_METRICS.sh"
                elif [ "$SEQ_TYPE" = "SE" ]
                    then
                        echo "Collecting single-end alignment metrics"
                        ALIGNMENT_METRICS_SCRIPT="9_SE_ALIGNMENT_METRICS.sh"
                else
                    echo -e "${red}ERROR - SEQ_TYPE parameter not recognised:${NC} "$SEQ_TYPE"
                    "
                    exit 1
                fi

                # sh "$PIPELINE_SCRIPTS_DIR"/9_PE_ALIGNMENT_METRICS.sh
                # Usage: 9_PE_ALIGNMENT_METRICS.sh <ANALYSIS_DIR> <QC_DIR> <SAMPLE_DIR> <SAMPLE_NAME> <PICARD_MEM> <REF_FILE> <PICARD_SIF> <SAMTOOLS_SIF>

                sbatch -W \
                        --account="$THIS_USER_ACCOUNT" \
                        --job-name="$JOB_NAME" \
                        --output="$STAGE_OUTPUT".%j.%N.out \
                        --error="$STAGE_ERROR".%j.%N.err \
                        --partition=defq \
                        --time=10:00:00 \
                        --nodes=1 \
                        --ntasks=1 \
                        --cpus-per-task=4 \
                        --mem-per-cpu=16G \
                        --wrap="\
                                sh "$PIPELINE_SCRIPTS_DIR"/"$ALIGNMENT_METRICS_SCRIPT" \
                                "$THIS_ANALYSIS_DIR" \
                                "$QC_DIR_NAME" \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ \
                                "$THIS_SAMPLE_NAME" \
                                "$PICARD_MEM_9" \
                                "$REF_FILE" \
                                "$PICARD_SIF" \
                                "$SAMTOOLS_SIF"\
                                "

                # Catch output status
                OUTPUT_STATUS=$?

                # Get Job ID
                getJobID "$JOB_NAME"

                # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                # Check for exit code error: <function> <$JOB_ID>
                slurmExitCodeCheck "$JOB_ID"
                
                # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"

                # check output status: <function> <stage label>
                checkOutputStatus "$STAGE_NAME"
       
                echo -e "${green}"$JOB_NAME" complete.${NC}
                "
                # update the stage
                START_AT_STAGE=$(($START_AT_STAGE+1))
                STAGE_NAME=""
                STAGE_OUTPUT=""
                STAGE_ERROR=""
        fi

        # check if we need to exit
        if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
        then
                echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
                exit 0
        fi


# 10) RSEQC
        # run the stage
        STAGE_NAME="10-RSEQC"
        if [ $START_AT_STAGE -eq 10 ]
        then
                JOB_NAME=stage-"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                echo -e "${blue}Running "$JOB_NAME"${NC}
                "
                
                # sh 10_RSEQC.sh
                # Usage: 10_RSEQC.sh <SEQ_TYPE> <ANALYSIS_DIR> <QC_DIR> <SAMPLE_DIR> <SAMPLE_NAME> <REFSEQ_BED> <HOUSEKEEPING_BED> <RSEQC_SIF> <RSEQC_REFS>

                sbatch -W \
                        --account="$THIS_USER_ACCOUNT" \
                        --job-name="$JOB_NAME" \
                        --output="$STAGE_OUTPUT".%j.%N.out \
                        --error="$STAGE_ERROR".%j.%N.err \
                        --partition=defq \
                        --time=10:00:00 \
                        --nodes=1 \
                        --ntasks=1 \
                        --cpus-per-task=8 \
                        --mem-per-cpu=16G \
                        --wrap="\
                                sh "$PIPELINE_SCRIPTS_DIR"/10_RSEQC.sh \
                                "$SEQ_TYPE" \
                                "$THIS_ANALYSIS_DIR" \
                                "$QC_DIR_NAME" \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ \
                                Alignments \
                                "$THIS_SAMPLE_NAME".mapped.rgid.filtered.sorted.bam \
                                "$THIS_SAMPLE_NAME" \
                                "$REFSEQ_BED" \
                                "$HOUSEKEEPING_BED" \
                                "$RSEQC_SIF" \
                                "$RSEQC_REFS"\
                                "

                # Catch output status
                OUTPUT_STATUS=$?

                # Get Job ID
                getJobID "$JOB_NAME"

                # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                # Check for exit code error: <function> <$JOB_ID>
                slurmExitCodeCheck "$JOB_ID"
                
                # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"

                # check output status: <function> <stage label>
                checkOutputStatus "$STAGE_NAME"
       
                echo -e "${green}"$JOB_NAME" complete.${NC}
                "
                # update the stage
                START_AT_STAGE=$(($START_AT_STAGE+1))
                STAGE_NAME=""
                STAGE_OUTPUT=""
                STAGE_ERROR=""
        fi

        # check if we need to exit
        if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
        then
                echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
                exit 0
        fi


# 11) HTSEQ COUNT
        # run the stage
        STAGE_NAME="11-HTSEQ-COUNT"
        if [ $START_AT_STAGE -eq 11 ]
        then
                JOB_NAME=stage-"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"
                STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"
                echo -e "${blue}Running "$JOB_NAME"${NC}
                "

                # sh 11_HTSEQ_COUNT.sh
                # Usage: 11_HTSEQ_COUNT.sh <ANALYSIS_DIR> <SAMPLE_DIR> <SAMPLE_NAME> <STRAND_TYPE> <SPIKE_IN> <SORT_ORDER> <GTF> <DM_GTF> <ERCC_GTF> <FEATURETYPE> <IDATTR> <MINAQUAL> <MODE> <COUNTS_DIRNAME> <SEQ_TYPE> <HTSEQ_SIF> <REFS_DIR>

                sbatch -W \
                        --account="$THIS_USER_ACCOUNT" \
                        --job-name="$JOB_NAME" \
                        --output="$STAGE_OUTPUT".%j.%N.out \
                        --error="$STAGE_ERROR".%j.%N.err \
                        --partition=defq \
                        --time=10:00:00 \
                        --nodes=1 \
                        --ntasks=1 \
                        --cpus-per-task=8 \
                        --mem-per-cpu=8G \
                        --wrap="\
                                sh "$PIPELINE_SCRIPTS_DIR"/11_HTSEQ_COUNT.sh \
                                "$THIS_ANALYSIS_DIR" \
                                "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ \
                                "$THIS_SAMPLE_NAME" \
                                "$STRAND_TYPE" \
                                "$SPIKE_IN" \
                                "$ORDER" \
                                "$GTF" \
                                "$DM_GTF" \
                                "$ERCC_GTF" \
                                "$FEATURETYPE" \
                                "$IDATTR" \
                                "$HT_MINAQUAL" \
                                "$MODE" \
                                "$COUNTS_DIRNAME" \
                                "$SEQ_TYPE" \
                                "$HTSEQ_SIF" \
                                "$REFS_DIR"\
                                "

                # Catch output status
                OUTPUT_STATUS=$?

                # Get Job ID
                getJobID "$JOB_NAME"

                # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
                slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

                # Check for exit code error: <function> <$JOB_ID>
                slurmExitCodeCheck "$JOB_ID"
                
                # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
                copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"

                # check output status: <function> <stage label>
                checkOutputStatus "$STAGE_NAME"
       
                echo -e "${green}"$JOB_NAME" complete.${NC}
                "
                # update the stage
                START_AT_STAGE=$(($START_AT_STAGE+1))
                STAGE_NAME=""
                STAGE_OUTPUT=""
                STAGE_ERROR=""
        fi

        # check if we need to exit
        if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
        then
                echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
                exit 0
        fi

# NOT IN USE for RNAseq
# # 12) HOMER Make Tag directory
#         # run the stage
#         STAGE_NAME="12-HOMER_MAKE_TAG_DIRECTORY"
#         if [ $START_AT_STAGE -eq 12 ]

#         then
#                 STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$PIPELINE_TYPE"-stage"$STAGE_NAME"-"$THIS_SAMPLE_NAME"
#                 STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$PIPELINE_TYPE"-stage"$STAGE_NAME"-"$THIS_SAMPLE_NAME"
#                 echo -e "${blue}Starting stage "$STAGE_NAME" for "$THIS_SAMPLE_NAME"${NC}
#                 "

#                 # sh "$PIPELINE_SCRIPTS_DIR"/12_HOMER_BAM_makeTagDirectory.sh
#                 # Usage: 12_HOMER_BAM_makeTagDirectory.sh <SAMPLE_DIR> <HOMER_DIRNAME> <BAM_IN_FILENAME> <SAMPLE_NAME> <STRAND_TYPE (strand-specific-fwd/strand-specific-rev/unstranded)> <SEQ_TYPE (PE/SE)> <REF> <FRAG_LENGTH> <TAGS_PER>

#                 sbatch -W \
#                         --account="$THIS_USER_ACCOUNT" \
#                         --job-name=stage"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME" \
#                         --output="$STAGE_OUTPUT".%j.%N.out \
#                         --error="$STAGE_ERROR".%j.%N.err \
#                         --time=24:00:00 \
#                         --mem=108G \
#                         --nodes=1 \
#                         --cpus-per-task=8 \
#                         --ntasks=1 \
#                         --wrap="\
#                                   sh "$PIPELINE_SCRIPTS_DIR"/12_HOMER_MAKE_TAG_DIR.sh "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ "$HOMER_DIRNAME" "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/Alignments/"$THIS_SAMPLE_NAME".mapped.rgid.sprted.filtereddups_mark.bam "$THIS_SAMPLE_NAME" "$STRAND_TYPE" "$SEQ_TYPE" "$REF" "$FRAG_LENGTH_s11" "$TAGS_PER_s11""

#                 # Catch output status
#                 OUTPUT_STATUS=$?

#                 # Get Job ID
#                 getJobID "$JOB_NAME"

#                 # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
#                 slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

#                 # Check for exit code error: <function> <$JOB_ID>
#                 slurmExitCodeCheck "$JOB_ID"
                
#                 # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
#                 copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"

#                 # check output status: <function> <stage label>
#                 checkOutputStatus "$STAGE_NAME"
       
#                 echo -e "${green}"$JOB_NAME" complete.${NC}
#                 "
#                 # update the stage
#                 START_AT_STAGE=$(($START_AT_STAGE+1))
#                 STAGE_NAME=""
#                 STAGE_OUTPUT=""
#                 STAGE_ERROR=""
#         fi

#         # check if we need to exit
#         if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
#         then
#                 echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
#                 exit 0
#         fi

# # 13) HOMER make UCSCfile
#         # run the stage
#         STAGE_NAME="13-HOMER_MAKE_UCSC_FILE"
#         if [ $START_AT_STAGE -eq 13 ]

#         then
#                 STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$PIPELINE_TYPE"-stage"$STAGE_NAME"-"$THIS_SAMPLE_NAME"
#                 STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$PIPELINE_TYPE"-stage"$STAGE_NAME"-"$THIS_SAMPLE_NAME"
#                 echo -e "${blue}Starting stage "$STAGE_NAME" for "$THIS_SAMPLE_NAME"${NC}
#                 "

#                 # sh "$PIPELINE_SCRIPTS_DIR"/13_HOMER_MAKE_UCSC_FILE.sh
#                 # Usage: 13_HOMER_MAKE_UCSC_FILE.sh <SAMPLE_DIR> <HOMER_TAG_DIRNAME> <SAMPLE_NAME> <TRACKS_DIRNAME> <STRAND_TYPE (strand-specific-fwd/strand-specific-rev/unstranded)> <REF> <FRAG_LENGTH> <TAGS_PER> <NORM> <RES>

#                 sbatch -W \
#                         --account="$THIS_USER_ACCOUNT" \
#                         --job-name=stage"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME" \
#                         --output="$STAGE_OUTPUT".%j.%N.out \
#                         --error="$STAGE_ERROR".%j.%N.err \
#                         --time=24:00:00 \
#                         --mem=108G \
#                         --nodes=1 \
#                         --cpus-per-task=8 \
#                         --ntasks=1 \
#                         --wrap="\
#                                 sh "$PIPELINE_SCRIPTS_DIR"/13_HOMER_MAKE_UCSC_FILE.sh "$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/ "$HOMER_TAG_DIRNAME" "$THIS_SAMPLE_NAME" "$TRACKS_DIRNAME" "$STRAND_TYPE" "$REF" "$FRAG_LENGTH_s12" "$TAGS_PER_s12" "$NORM_s12" "$RES_s12""

#                 # Catch output status
#                 OUTPUT_STATUS=$?

#                 # Get Job ID
#                 getJobID "$JOB_NAME"

#                 # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
#                 slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"

#                 # Check for exit code error: <function> <$JOB_ID>
#                 slurmExitCodeCheck "$JOB_ID"
                
#                 # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
#                 copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"

#                 # check output status: <function> <stage label>
#                 checkOutputStatus "$STAGE_NAME"
       
#                 echo -e "${green}"$JOB_NAME" complete.${NC}
#                 "
#                 # update the stage
#                 START_AT_STAGE=$(($START_AT_STAGE+1))
#                 STAGE_NAME=""
#                 STAGE_OUTPUT=""
#                 STAGE_ERROR=""
#         fi

#         # check if we need to exit
#         if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
#         then
#                 echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
#                 exit 0
#         fi

##########################################
# # NN) Stage Template    ## (some stages might require detection of SE/PE etc eg see FASTQC stages)  ## UNCOMMENT WHOLE SECTION AND REMOVE ALL COMMENTS WITH "##" BEFORE USING
#         # run the stage
#         STAGE_NAME="NN-STAGE-NAME"      ## Change NN to suit
#         if [ $START_AT_STAGE -eq NN ]   ## Change NN to suit
#         then
#                 JOB_NAME=stage-"$STAGE_NAME"-"$PIPELINE_TYPE"-"$THIS_SAMPLE_NAME"    ## This is used by several other parts -can add extra labels in here eg "_R1"
#                 STAGE_OUTPUT="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"       ## Do not change: this is used for creation of log file and cleanup
#                 STAGE_ERROR="$THIS_ANALYSIS_DIR"/Sample_"$THIS_SAMPLE_NAME"/std_out_err/"$JOB_NAME"       ## Do not change: this is used for creation of log file and cleanup
#                 echo -e "${blue}Running "$JOB_NAME"${NC}
#                 "   ## Do not change: This echo line is used for finding job ID
                
#                 # sh NN_STAGE_NAME.sh
#                 # Usage: NN_STAGE_NAME.sh <ARG1> <ARG2> <ARG3> <ARG4> <etc...>  ## ARGS will use global and/or stage-specific variables entered in section at top of this pipeline script
#                                                                                 ## Any variables specific to this stage should be added to editable variables section at top to facilitate easy customization

#                 sbatch -W \     ## -W forces sbatch to wait for command inside --wrap to finish
#                         --account="$THIS_USER_ACCOUNT" \       ## Do not change
#                         --job-name="$JOB_NAME" \       ## Do not change
#                         --output="$STAGE_OUTPUT".%j.%N.out \       ## Do not change: this is used for creation of log file and cleanup
#                         --error="$STAGE_ERROR".%j.%N.err \       ## Do not change: this is used for creation of log file and cleanup
#                         --partition=defq \  ## Can use bigmem for some stages
#                         --time=10:00:00 \    ## Alter to suit
#                         --mem=32GB \    ## Alter to suit
#                         --nodes=1 \       ## Do not change
#                         --cpus-per-task=8 \     ## Alter to suit
#                         --ntasks=1 \       ## Do not change
#                         --wrap="\
#                                 sh "$PIPELINE_SCRIPTS_DIR"/NN_STAGE_NAME.sh <ARG1> <ARG2> <ARG3> <ARG4> <etc...>"   ## --wrap sends contents as shell script to slurm. This is where the actual stage script is called and ARGS are passed to be used as variables. Can be broken onto separate lines with "\" for ease of reading but be careful of formatting.

#                 ## Do not remove the following lines as they are used to catch errors and compile log files
#                 # Catch output status
#                 OUTPUT_STATUS=$?       ## Do not change

#                 # Get Job ID
#                 getJobID "$JOB_NAME"       ## Do not change

#                 # Check for SLURM error: <function> <$STAGE_ERROR> <$JOB_ID>
#                 slurmErrorCheck "$STAGE_ERROR" "$JOB_ID"       ## Do not change

#                 # Check for exit code error: <function> <$JOB_ID>
#                 slurmExitCodeCheck "$JOB_ID"       ## Do not change
                
#                 # Copy standard output and standard error to logfile: <function> <Title> <$STAGE_ERROR> <$JOB_ID>
#                 copyOutErrLog "$STAGE_NAME" "$STAGE_ERROR" "$JOB_ID"       ## First ARG <Title> can have extra text added eg "_R1" if needed; otherwise leave as $STAGE_NAME

#                 # check output status: <function> <stage label>       ## First ARG <stage label> can have extra text added eg "_R1" if needed; otherwise leave as $STAGE_NAME
#                 checkOutputStatus "$STAGE_NAME"
       
#                 echo -e "${green}"$JOB_NAME" complete.${NC}\n"       ## Do not change
#                 # update the stage
#                 START_AT_STAGE=$(($START_AT_STAGE+1))       ## Do not change
#                 STAGE_NAME=""       ## Do not change
#                 STAGE_OUTPUT=""       ## Do not change
#                 STAGE_ERROR=""       ## Do not change
#         fi

#         # check if we need to exit       ## Do not change
#         if [ $END_AT_STAGE -eq $(($START_AT_STAGE-1)) ]
#         then
#                 echo -e "------\n${green}"$PIPELINE_TITLE" run ending at stage "$END_AT_STAGE"${NC}\n------\n"
#                 exit 0
#         fi
##########################################




# check output status
if [ $? -ne 0 ]
then
        echo -e "     ----------\n${green}All "$PIPELINE_TYPE" stages run and completed for sample: "$THIS_SAMPLE_NAME"\n----------"
fi

exit







