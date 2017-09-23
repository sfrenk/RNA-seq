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

#       default directory:
dir="."

usage="
    Assemble de-novo transcriptomes for samples with the stringtie pipeline. NOTE: This pipeline requires mapped bam files. If starting with raw data (eg. .fastq files) run hisat2.sh first.

    ARGUMENTS
        -d/--dir
        Master directory containing the following directories:
            ./bam/
            ./stringtie/

        (default: current directory)

    "

# Parse command line parameters

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -d|--dir)
        dir="$2"
        shift
        ;;
        -h|--help)
        echo "$usage"
        ;;
    esac
shift
done

###############################################################################
###############################################################################

if [[ -f ./stringtie/stringtie_merged.gtf ]]; then
    rm ./stringtie/stringtie_merged.gtf
fi

ls ${dir}/stringtie/*.gtf > gtf_files.txt

stringtie --merge -p $SLURM_NTASKS -G $gtf -o ./stringtie/stringtie_merged.gtf gtf_files.txt

for file in ${dir}/bam/*.bam; do
	base=$(basename $file .bam)

    stringtie $file -p $SLURM_NTASKS -G ./stringtie/stringtie_merged.gtf -e -B -o ./stringtie/${base}/transcripts.gtf

done

