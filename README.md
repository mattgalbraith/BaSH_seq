# BaSH_seq Pipeline (**B**ash **a**nd **S**lurm on HPC)  

Developed and maintained by Matthew Galbraith  

This pipeline processes high-throughput sequencing data (PE/SE) through QC,trim/filter,alignment,counting etc via sequential stages, with individual samples run in parallel via submission to a Slurm queue.  
This specific version was customized to run on HDC Eureka instances (running CentOS 7) on Google Compute Engine and uses a custom Conda environment to supply most of the tools/dependencies required by each stage, with only a few problemmatic tools running from Singularity containers.


```
BaSH_seq/
├── CutRun
│   └── stageScripts  # Cut&Run-specific stage scripts + symbolic links to common stage scripts
├── RNAseq
│   └── stageScripts  # RNAseq-specific stage scripts + symbolic links to common stage scripts
├── misc
├── multiqc
├── other_scripts
├── pipeline_templates
└── stageScripts                # common stage scripts live here
```


## Dependencies:
Bash
Slurm
Tools required by each stage (in PATH / Conda env / Containers); Conda / Docker/Singularity
References and related files
MultiQC script (need to update and integrate)
R (if integrating RPKM script)

## Installation and setup 
need to download required refs etc 
need to modify paths to stage scripts, references, etc
fastq_screen.conf
For Conda version: need to add env.yaml
For container version: modify stage scripts to use correct call to each tool or container



## Pre-run steps?
FASTQ merging?


## Steps to run pipeline:  
1) In `Project/` : create top-level working dirs: `Project/raw_date` and `Project/analysis_date`  
2) In `Project/analysis_date/` create sample_locations.txt (field 1 = SAMPLE_NAME; field 2 = path/to/raw_fastq_file.gz; field 3 = read1 / read2 labels)  
    See also: SampleInfo.xlsx template
3) Edit variables in top section (SEE XXX FOR MORE DETAILS eg Human vs Mouse & strandedness) of pipeline script template and save a project-specific version in `Project/analysis_date/scripts`(specify path to this version as indicated below)  
4) From `Project/analysis_date/` run analysis_setup.sh as follows:  
```sh path/to/analysis_setup.sh <FULL/PATH/TO/PIPELINE_SCRIPT> <START_AT_STAGE> <END_AT_STAGE>```  
	This will create and populate `Project/analysis_date/Sample_*` directories with symbolic links to FASTQ files in `Project/raw_date` as well as write out 'submit' scripts to `Project/analysis_date/scripts`  
	Usually will specify subsets of pipeline stages to allow for QC checks/troubleshooting eg stages 1-3, 4-4, 4-11  
5) Then from `Project/analysis_date/scripts` run submit scripts as follows (command can be copied from within submit script):  
```sbatch submitAll_START_AT_STAGE-END_AT_STAGE.sh```  
    If using conda to magae required tools, will need to first activate appropriate env  
6) Monitor progress using `squeue` (see watch alias)


## Post-run steps and troubleshooting  
1) check pipeline logs
2) (if needed) check Sample/Stage logs
3) Run MultiQC and inspect report
4) (RNAseq only) Gather counts files 
5) Run temp FASTQ and BAM file cleanup
6) (if needed) sync back to main data store