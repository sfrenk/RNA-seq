#!/usr/bin/env bash

###############################################################################
# Basic gene expression analysis pipeline for mapping single or paired end fastq reads to mRNA or transposon consensus sequences.
###############################################################################

# This script produces the following output files:
# Indexed bam file for each sample
# Count table for all samples combined


# Before running the script, make sure modules are loaded:
#       bbmap
#       bowtie2
#       samtools
#       python (default version)

# The location of the following files may have to be modified in this script:
# bowtie2 index files
# merge_counts.py
# illumina adaptor sequences

# NOTE: Any read processing beyond adapter trimming is not covered in this pipeline and must be performed beforehand.

###############################################################################
###############################################################################

usage="
    USAGE
       step1:   load the following modules: bbmap bowtie2 samtools python (default version)
       step2:   bash bowtie2_pipeline.sh [options]  

    ARGUMENTS
        -f/--fastq
        Directory containing fastq files (these files should be gzipped)
        
        -r/--reference
        Mapping reference (mrna or transposons)

        -p/--paired
        Use this option if fastq files contain paired-end reads
    "
# Set default reference

REF="mrna"

# Parse command line parameters

if [ -z "$1" ]; then
        echo "$usage"
        exit
fi

while [[ $# > 0 ]]
do
        key="$1"
        case $key in
                -f|--fastq)
                FASTQ="$2"
                shift
                ;;
                -r|--ref)
                REF="$2"
                shift
                ;;
                -p|--paired)
                PAIRED="true"
                ;;
        esac
shift
done

# Remove trailing "/" from fastq directory if present

if [[ ${FASTQ:(-1)} == "/" ]]; then
        FASTQ=${FASTQ::${#FASTQ}-1}
fi


# Select bowtie2 index based on reference

if [[ $REF == "mrna" ]]; then
        INDEX="/nas02/home/s/f/sfrenk/proj/seq/WS251/mrna/bowtie2/mrna"
elif [[ $REF == "transposons" ]]; then
        INDEX="/proj/ahmedlab/steve/seq/transposons/bowtie2/transposon"
else
        echo "ERROR: Invalid reference argument"
fi

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
if [ ! -d "fastq_trimmed" ]; then
        mkdir fastq_trimmed
fi

# Print out loaded modules to keep a record of which software versions were used in this run

modules=$(/nas02/apps/Modules/bin/modulecmd tcsh list 2>&1)
echo "$modules"

# Start pipeline

echo $(date +"%m-%d-%Y_%H:%M")" Starting pipeline..."

for file in ${FASTQ}/*.fastq.gz
do
        if [[ $PAIRED == "true" ]]; then
            
            # Paired end

            if [[ ${file:(-11)} == "r1.fastq.gz" ]]; then

                # Process and map r1 and r2 reads simultaniously

                FBASE=$(basename $file .fastq.gz)
                BASE=${FBASE%_r1}

                # Remove adaptor sequence

                echo $(date +"%m-%d-%Y_%H:%M")"Processing ${BASE}"
                bbduk.sh in1=${FASTQ}/${BASE}_r1.fastq.gz in2=${FASTQ}/${BASE}_r2.fastq.gz out1=./fastq_trimmed/${BASE}_r1.fastq out2=./fastq_trimmed/${BASE}_r2.fastq ktrim=r k=23 hdist=1 ref=/nas02/apps/bbmap-34.83/bbmap/resources/truseq_rna.fa.gz tpe tbo
                
                echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${BASE} with Bowtie2..."

                # Map reads using Bowtie 2

                bowtie2 -q -S ./bowtie2_out/${BASE}.sam -p 4 --no-mixed -x ${INDEX} -1 ./fastq_trimmed/${BASE}_r1.fastq -2 ./fastq_trimmed/${BASE}_r2.fastq
            else

                # Avoid double mapping by skipping the r2 read file
                
                SKIPFILE="true"
            fi
        else
            
            # Single end

            BASE=$(basename $file .fastq.gz)

            echo $(date +"%m-%d-%Y_%H:%M")"Processing ${BASE}"

            # Remove adaptor sequence

            echo $(date +"%m-%d-%Y_%H:%M")" Removing adaptor sequences..."
            bbduk.sh in=$file out=./fastq_trimmed/${BASE}.fastq ktrim=r k=23 hdist=1 ref=/nas02/apps/bbmap-34.83/bbmap/resources/truseq_rna.fa.gz

            # Map reads using Bowtie 2

            echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${BASE} with Bowtie2..."

            bowtie2 -q -S ./bowtie2_out/${BASE}.sam -p 4 -x ${INDEX} -U ./fastq_trimmed/${BASE}.fastq
        fi     
        
        if [[ $SKIPFILE != "true" ]]; then
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

        else
            SKIPFILE="false"
        fi

done

# Create count table

echo $(date +"%m-%d-%Y_%H:%M")"Merging count files into count table"

python /proj/ahmedlab/steve/seq/util/merge_counts.py ./count
