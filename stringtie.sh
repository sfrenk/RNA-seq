#!/usr/bin/env bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -t 1-0

module add stringtie
module list

###############################################################################
# Assemble de-novo transcriptomes for samples with the cufflinks pipeline.
###############################################################################

# VARIABLES
#       gtf annotation file for genes
gtf="/nas02/home/s/f/sfrenk/proj/seq/WS251/genes.gtf"

usage="
    Assemble de-novo transcriptomes for samples with the stringtie pipeline. NOTE: This pipeline requires mapped bam files. If starting with raw data (eg. .fastq files) run hisat2.sh first.

    ARGUMENTS
        -d/--dir
        Master directory containing the following directories:
            ./bam/
            ./stringtie/

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

shopt -s nullglob

gtf_files=(./stringtie/*.gtf)

stringtie --merge -p $SLURM_NTASKS -G $gtf -o ./stringtie/stringtie_merged.gtf ${gtf_files[@]}

for file in ${dir}/*.bam; do
	base=$(basename $file .bam)

    stringtie $file -p $SLURM_NTASKS -G ./stringtie/stringtie_merged.gtf -e -b ./stringtie/${base}.ctab


done

