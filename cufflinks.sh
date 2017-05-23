#!/usr/bin/env bash
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -t 7-0

module add samtools cufflinks

###############################################################################
# Assemble de-novo transcriptomes for samples with the cufflinks pipeline.
###############################################################################



usage="
    Assemble de-novo transcriptomes for samples with the cufflinks pipeline. NOTE: This pipeline requires mapped bam files. If starting with raw data (eg. .fastq files) run tophat2.sh first.

    USAGE
       step1:   load the following modules: trim_galore hisat2 samtools subread
       step2:   bash hisat2_genome.sh [options]  

    ARGUMENTS
        -d/--dir
        Directory containing bam files.

    "

# Parse command line parameters

if [ -z "$1" ]; then
    echo "$usage"
    exit
fi


while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -d|--dir)
        dir="$2"
        shift
        ;;
    esac
shift
done

# Remove trailing "/" from input directory if present

if [[ ${dir:(-1)} == "/" ]]; then
    dir=${dir::${#dir}-1}
fi


###############################################################################
###############################################################################


# Make directories

if [ ! -d "sam" ]; then
    mkdir sam
fi

for file in ${dir}/*.bam; do
	base=$(basename $file .bam)
	
	samtools view -f $file > ${base}.sam

	cufflinks -p 4 ${base}.sam

	cuffmerge -g genes.gtf -s genome -p 4 assemblies.tx
