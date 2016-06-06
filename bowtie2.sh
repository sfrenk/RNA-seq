#!/usr/bin/env bash

###############################################################################
# Basic gene expression analysis pipeline for mapping single or paired end fastq reads to mRNA or transposon consensus sequences.
###############################################################################

# This script produces the following output files:
# Indexed, sorted bam file for each sample
# Count table for all samples combined


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
        Directory containing fastq files (these files should be gzipped)
        
        -r/--reference
        Mapping reference (mrna or transposons)

        -p/--paired
        Use this option if fastq files contain paired-end reads. NOTE: if paired, each pair must consist of two files with the basename ending in '_r1' or '_r2' depending on respective orientation.
    "
# Set default reference

REF="mrna"
PAIRED=false

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
        DIR="$2"
        shift
        ;;
        -r|--ref)
        REF="$2"
        shift
        ;;
        -p|--paired)
        PAIRED=true
        ;;
    esac
shift
done

# Remove trailing "/" from fastq directory if present

if [[ ${DIR:(-1)} == "/" ]]; then
    DIR=${DIR::${#DIR}-1}
fi


# Select bowtie2 index based on reference

if [[ $REF == "mrna" ]]; then
    INDEX="/nas02/home/s/f/sfrenk/proj/seq/WS251/mrna/bowtie2/mrna"
elif [[ $REF == "transposons" ]]; then
    INDEX="/proj/ahmedlab/steve/seq/transposons/bowtie2/transposon"
else
    echo "ERROR: Invalid reference argument"
fi

# Print out loaded modules to keep a record of which software versions were used in this run

echo "$modules"

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

for file in ${DIR}/*.fastq.gz; do

    SKIPFILE=false
        
    if [[ $PAIRED = true ]]; then
            
        # Paired end

        if [[ ${file:(-11)} == "r1.fastq.gz" ]]; then

            # Process and map r1 and r2 reads simultaniously

            FBASE=$(basename $file .fastq.gz)
            BASE=${FBASE%_r1}
                
            echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${BASE} with Bowtie2..."

            # Map reads using Bowtie 2

            bowtie2 -q -S ./bowtie2_out/${BASE}.sam -p 4 --no-mixed -x ${INDEX} -1 ${DIR}/${BASE}_r1.fastq.gz -2 ${DIR}/${BASE}_r2.fastq.gz
        else

            # Avoid double mapping by skipping the r2 read file and proceding to the next step 
                
            SKIPFILE=true
        fi
    else
            
        # Single end

        BASE=$(basename $file .fastq.gz)

        # Map reads using Bowtie 2

        echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${BASE} with Bowtie2..."

        bowtie2 -q -S ./bowtie2_out/${BASE}.sam -p 4 -x ${INDEX} -U $file
    fi     
        
    if [[ $SKIPFILE = false ]]; then
        echo $(date +"%m-%d-%Y_%H:%M")" Mapped ${BASE}"

        # Convert to bam then sort

        echo $(date +"%m-%d-%Y_%H:%M")" Converting and sorting ${BASE}..."

        samtools view -bS ./bowtie2_out/${BASE}.sam > ./bam/${BASE}.bam

        samtools sort -o ./bam/${BASE}_sorted.bam ./bam/${BASE}.bam 

        # Need to index the sorted bam files for visualization

        echo $(date +"%m-%d-%Y_%H:%M")" Indexing ${BASE}..."

        samtools index ./bam/${BASE}_sorted.bam
        
        # Count reads

        echo $(date +"%m-%d-%Y_%H:%M")" Counting reads"
        grep -v '^@' bowtie2_out/${BASE}.sam | cut -f 3 | sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' > ./count/${BASE}_counts.txt

        # Remove sam file

        rm ./bowtie2_out/${BASE}.sam
    fi

done

# Create count table

echo $(date +"%m-%d-%Y_%H:%M")"Merging count files into count table"

python /proj/ahmedlab/steve/seq/util/merge_counts.py ./count
