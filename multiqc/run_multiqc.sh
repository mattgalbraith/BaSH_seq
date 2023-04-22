#!/bin/bash

#run with ./run_multiqc.sh /Users/kinninko/Projects/HCT116_Hypoxia_spikein/analysis_05_06_2019 brief

SCRIPT_TITLE=run_multiqc.sh
SCRIPT_VERSION=0.1
# DATE: 11-08-2019
# AUTHOR: Kohl Kinning
# SUMMARY: 
# This script can be run from anywhere and will create the directories:
# Project/analysis_mm_dd_yyyy//QC/multiqc/ and
# Project/analysis_mm_dd_yyyy//QC/multiqc/summary_counts/
# and will generate:
# Project/analysis_mm_dd_yyyy/QC/multiqc/multiqc_report_[brief|verbose].html
# and some tables in Project/analysis_mm_dd_yyyy/QC/multiqc/summary_counts/
# 
# This script will scan the analysis folder for log files depending on verbosity requested and use the standard config files stored in $SHARED


##### ARGUMENTS PASPED FROM COMMAND LINE #####

# variables from command line:
WD=${1}
VERBOSITY=${2}

EXPECTED_ARGS=2
# check if correct number of arguments are supplied from command line
if [ $# -ne $EXPECTED_ARGS ]
then
        echo -ne "${red}ERROR - expecting "$EXPECTED_ARGS" ARGS but "$#" provided:${NC}
        "$@""
        echo "Usage: sh <PATH/TO/run_multiqc.sh> <PATH/TO/ANALYSIS/FOLDER> <VERBOSITY: brief/verbose>"
        echo -e "EXAMPLE: sh /gpfs/jespinosa/shared/PipeLinesScripts/run_multiqc.sh /gpfs/jespinosa/shared/Projects2/RNAseq/HCT116_Hypoxia_spikein/analysis_05_06_2019 brief"
        exit 1
fi


# CHANGEPATHS
if [ $VERBOSITY = "brief" ]; then
  CONFIG_FILE="/Users/kohl/Projects/pipeline_scripts/multiqc/multiqc_config_brief.yaml"
  elif [ $VERBOSITY = "verbose" ]; then
  CONFIG_FILE="/Users/kohl/Projects/pipeline_scripts/multiqc/multiqc_config_verbose.yaml"
fi

# Check for sample_locations.txt
if [ ! -e $WD/sample_locations.txt ]
	then
		echo -e "${red}ERROR - sample_locations.txt not found
		"
		exit 1
fi

#helper function for joining with a delimiter
function join_by { local IFS="$1"; shift; echo "$*"; }

#create folder (if it DNE) and initialize table with header
mkdir -p $WD/QC/multiqc/
mkdir -p $WD/QC/multiqc/summary_counts/
echo "Sample_name,Reads_pre_trim,Reads_post_trim,Aligned_reads_pre_filter,Aligned_reads_post_filter,Duplicate_reads,HTSeq_counts" > $WD/QC/multiqc/summary_counts/final_read_summary.csv

#grep sample name from file
for i in `grep -e read1 $WD/sample_locations.txt | awk '{print $1}'`;
do
    if [[ ! -f $WD/QC/PE_metrics/"$i"_PE_alignmentMetricsSummary.txt ]] ; then
    echo 'PE_alignmentMetricsSummary.txt does not exist. You will not have any custom plots or tables.'
    #exit 1
    fi

    counts[0]=$i
    #read the stats file 
    counts[1]=$(cat $WD/QC/PE_metrics/"$i"_PE_alignmentMetricsSummary.txt | grep 'Total Sequences' | head -1 | awk '{print $4}');
    counts[2]=$(cat $WD/QC/PE_metrics/"$i"_PE_alignmentMetricsSummary.txt | grep 'Total Sequences' | tail -1 | awk '{print $4}');
    #halve the count
    counts[3]=$(($(cat $WD/QC/PE_metrics/"$i"_PE_alignmentMetricsSummary.txt | grep "Aligned" | head -1 | awk '{print $2}') / 2));
    #halve the count
    counts[4]=$(($(cat $WD/QC/PE_metrics/"$i"_PE_alignmentMetricsSummary.txt | grep "Aligned" | tail -1 | awk '{print $2}') / 2));
    #halve the count
    counts[5]=$(($(cat $WD/QC/PE_metrics/"$i"_PE_alignmentMetricsSummary.txt | grep "Duplicates" | tail -1 | awk '{print $2}') / 2));
    counts[6]=$(awk '{sum+=$2}END{print sum}' $WD/Sample_"$i"/Counts/HTSeq/"$i"_HTSeq_counts.txt)
    
join_by , "${counts[@]}" >> $WD/QC/multiqc/summary_counts/final_read_summary.csv
done

#insert YAML header for multiQC bargraph
echo "\
# id: 'custom_read_summary_graph'
# section_name: 'Read summary (graph)'
# description: 'before and after trimming and filtering. Click the legend to hide/display categories.'
# format: 'csv'
# plot_type: 'bargraph'
# pconfig:
#    id: 'read_summary_graph'
#    ylab: 'Total reads'
#    stacking: NULL
#    tt_percentages: False
$(cat $WD/QC/multiqc/summary_counts/final_read_summary.csv)" > $WD/QC/multiqc/summary_counts/final_summary_graph_mqc.csv

#insert YAML header for multiQC bargraph, simpler counts
echo "\
# id: 'custom_feature_counts_graph'
# section_name: 'Feature counts (graph)'
# description: 'from HTSeq, post alignment.'
# format: 'csv'
# plot_type: 'bargraph'
# pconfig:
#    id: 'feature_counts_graph'
#    ylab: 'Total counts'
#    stacking: NULL
#    tt_percentages: False
#    cp_switch: False
$(cut -d, -f1,7 < $WD/QC/multiqc/summary_counts/final_read_summary.csv)" > $WD/QC/multiqc/summary_counts/final_counts_graph_mqc.csv

#insert YAML for multiQC table
echo "\
# id: 'custom_read_summary_table'
# section_name: 'Read summary (table)'
# description: 'before and after trimming and filtering.'
# format: 'csv'
# plot_type: 'table'
$(cat $WD/QC/multiqc/summary_counts/final_read_summary.csv)" > $WD/QC/multiqc/summary_counts/final_summary_table_mqc.csv

#GREEN='\033[0;32m'
#NC='\033[0m'
#status=$?
#[ $status -eq 0 ] && echo -e "Files \"final_read_summary.csv\", \"final_summary_table_mqc.csv\", \"final_summary_bargraph_mqc.csv\", and \"final_counts_bargraph_mqc.csv\" ${GREEN}successfully${NC} generated!"

multiqc -f -ip -c ${CONFIG_FILE} -o $WD/QC/multiqc -n multiqc_report_$VERBOSITY $WD

exit
