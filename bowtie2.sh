#!/usr/bin/env bash

###############################################################################
# Basic gene expression analysis pipeline for mapping single or paired end fastq reads to mRNA or transposon consensus sequences.
###############################################################################

# This script produces the following output files:
#       indexed, sorted bam file for each sample
#       Count table for all samples combined


# Before running the script, make sure modules are loaded:
#       bowtie2
#       samtools
#       python (default version)

# Check the location of the following files may have to be modified in this script:
#       bowtie2 index files
#       merge_counts.py

#       Set the command to display software versions used during the run
modules=$(/nas02/apps/Modules/bin/modulecmd tcsh list 2>&1)


# NOTE: Read processing is NOT covered in this pipeline and must be performed beforehand.

###############################################################################
###############################################################################

usage="
    USAGE
       step1:   load the following modules: bbmap bowtie2 samtools python (default version)
       step2:   bash bowtie2_pipeline.sh [options]  

    ARGUMENTS
        -d/--dir
        directory containing fastq files (these files should be gzipped)
        
        -r/--reference
        Mapping reference (mrna or transposons)

        -p/--paired
        Use this option if fastq files contain paired-end reads. NOTE: if paired, each pair must consist of two files with the basename ending in '_r1' or '_r2' depending on respective orientation.
    "
# Set default reference

ref="mrna"
paired=false

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
        -r|--ref)
        ref="$2"
        shift
        ;;
        -p|--paired)
        paired=true
        ;;
    esac
shift
done

# Remove trailing "/" from fastq directory if present

if [[ ${dir:(-1)} == "/" ]]; then
    dir=${dir::${#dir}-1}
fi


# Select bowtie2 index based on reference

if [[ $ref == "mrna" ]]; then
    index="/nas02/home/s/f/sfrenk/proj/seq/WS251/mrna/bowtie2/mrna"
elif [[ $ref == "transposons" ]]; then
    index="/proj/ahmedlab/steve/seq/transposons/bowtie2/transposon"
else
    echo "ERROR: Invalid reference argument"
fi

# Print out loaded modules to keep a record of which software versions were used in this run

echo "$modules"

# module test

req_modules=("trim_galore" "tophat" "samtools" "subread")

for i in ${req_modules[@]}; do
    if [[ $modules != *${i}* ]]; then
        echo "ERROR: Please load ${i}"
        exit 1
    fi
done

###############################################################################
###############################################################################

# Make directories

if [ ! -d "bowtie2_out" ]; then
    mkdir bowtie2_out
fi
if [ ! -d "count" ]; then
    mkdir count
fi
if [ ! -d "bam" ]; then
    mkdir bam
fi

# Start pipeline

echo $(date +"%m-%d-%Y_%H:%M")" Starting pipeline..."

for file in ${dir}/*.fastq.gz; do

    # The skipfile variable ensures that paired end reads don't get counted twice 

    skipfile=false
        
    if [[ $paired = true ]]; then
            
        # paired end

        if [[ ${file:(-11)} == "_1.fastq.gz" ]]; then

            # Process and map r1 and r2 reads simultaniously

            Fbase=$(basename $file .fastq.gz)
            base=${Fbase%_1}
                
            echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${base} with Bowtie2..."

            # Map reads using Bowtie 2

            bowtie2 -q -S ./bowtie2_out/${base}.sam -p 4 --no-mixed -x ${index} -1 ${dir}/${base}_1.fastq.gz -2 ${dir}/${base}_2.fastq.gz
        else

            # Avoid double mapping by skipping the r2 read file and proceding to the next step 
                
            skipfile=true
        fi
    else
            
        # Single end

        base=$(basename $file .fastq.gz)

        # Map reads using Bowtie 2

        echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${base} with Bowtie2..."

        bowtie2 -q -S ./bowtie2_out/${base}.sam -p 4 -x ${index} -U $file
    fi     
        
    if [[ $skipfile = false ]]; then
        echo $(date +"%m-%d-%Y_%H:%M")" Mapped ${base}"

        # Convert to bam then sort

        echo $(date +"%m-%d-%Y_%H:%M")" Converting and sorting ${base}..."

        samtools view -bS ./bowtie2_out/${base}.sam > ./bam/${base}.bam

        samtools sort -o ./bam/${base}_sorted.bam ./bam/${base}.bam 

        # Need to index the sorted bam files for visualization

        echo $(date +"%m-%d-%Y_%H:%M")" indexing ${base}..."

        samtools index ./bam/${base}_sorted.bam
        
        # Count reads

        echo $(date +"%m-%d-%Y_%H:%M")" Counting reads"
        # The 260 flag marks unmapped reads/secondary mappings
        awk -F$'\t' '$2 != "260" ' ./bowtie2_out/${base}.sam | cut -f 3 | sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' > ./count/${base}_counts.txt

        # Remove sam and unsorted bam files

        rm ./bowtie2_out/${base}.sam
        rm ./bam/${base}.bam
    fi

done

# Create count table

echo $(date +"%m-%d-%Y_%H:%M")"Merging count files into count table"

python /proj/ahmedlab/steve/seq/util/merge_counts.py ./count
