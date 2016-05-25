#!/usr/bin/env bash

###############################################################################
# Basic gene expression analysis pipeline for single end data using bowtie2
###############################################################################

###############################################################################
# BEFORE RUNNING SCRIPT DO THE FOLLOWING:
###############################################################################

# 1. Make sure modules are loaded:
#       bbmap
#       bowtie2
#       python (default version)

# The location of the following files may have to be modified in this script:
# bowtie index files
# merge_counts.py and small_rna_filter.py
# illumina adaptor sequences

# NOTE: Any read processing INCLUDING adapter trimming is not covered in this pipeline and must be performed beforehand.

###############################################################################
###############################################################################



usage="
    USAGE
       step1:   load the following modules: bbmap bowtie2 samtools python (default version)
       step2:   bash bowtie2_pipeline.sh [options]  

    ARGUMENTS
        -d/--dir
        Directory containing read files (can be .DIR .fasta or .txt (raw) format)

        -f 
        filter for 22G RNAs ('g') 21U RNAs ('u') or both ('gu' or 'ug')

        -c/--count
        Make a read count table
    "
# Set default parameters

FILTER=""
FILTER_G=false
FILTER_U=false
COUNT=false

# Parse command line parameters

if [ -z "$1" ]; then
        echo "$usage"
        exit
fi

while [[ $# > 0 ]]
do
        key="$1"
        case $key in
                -d|--DIR)
                DIR="$2"
                shift
                ;;
                -f|filter)
                FILTER="$2"
                shift
                ;;
                -c|--count)
                COUNT="true"
                ;;
        esac
shift
done

# Remove trailing "/" from DIR directory if present

if [[ ${DIR:(-1)} == "/" ]]; then
        DIR=${DIR::${#DIR}-1}
fi

# parse filter option
case $FILTER in
        "g")
        FILTER_G=true
        FILTER_U=false
        ;;
        "u")
        FILTER_G=false
        FILTER_U=true
        ;;
        "gu"|"ug")
        FILTER_G=true
        FILTER_U=true
        ;;
esac

###############################################################################
###############################################################################

# Make directories
if [ ! -d "filtered" ]; then
        mkdir filtered
fi

if [ ! -d "bowtie_out" ]; then
        mkdir bowtie_out
fi

if [ ! -d "bam" ]; then
        mkdir bam
fi

if [ ! -d "count" ] && [ "$COUNT" = true ]; then
        mkdir count
fi

# Print out loaded modules to keep a record of which software versions were used in this run

modules=$(/nas02/apps/Modules/bin/modulecmd tcsh list 2>&1)
echo "$modules"

# Start pipeline

echo $(date +"%m-%d-%Y_%H:%M")" Starting pipeline..."

for file in ${DIR}/*
do
        FBASE=$(basename $file .txt)
        BASE=${FBASE%.*}

        if [ "$FILTER_G" = true ]; then
                
                # Extract 22G RNAs

                echo $(date +"%m-%d-%Y_%H:%M")" Extracting 22G RNAs from $file ..."
                
                python /proj/ahmedlab/steve/seq/util/small_rna_filter.py -g -o ./filtered/${BASE}_22g.txt $file

                # Map reads using bowtie
                
                echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${BASE} 22G RNAs with Bowtie..."
                
                bowtie -p 4 -r -S -v 0 -a --best --strata -m 400 /proj/ahmedlab/steve/seq/WS251/genome/bowtie/genome ./filtered/${BASE}_22g.txt ./bowtie_out/${BASE}_22g.sam
        
                # Convert to bam then sort

                echo $(date +"%m-%d-%Y_%H:%M")" Converting and sorting ${BASE}..."

                samtools view -bS ./bowtie_out/${BASE}_22g.sam | samtools sort -o ./bam/${BASE}_22g.bam -

                # Need to index the sorted bam files for visualization

                echo $(date +"%m-%d-%Y_%H:%M")" Indexing ${BASE}..."

                samtools index ./bam/${BASE}_22g.bam
        fi

        if [ "$FILTER_U" = true ]; then
                
                # Extract 21U RNAs

                echo $(date +"%m-%d-%Y_%H:%M")" Extracting 21U RNAs from $file ..."
                
                python /proj/ahmedlab/steve/seq/util/small_rna_filter.py -u -o ./filtered/${BASE}_21u.txt $file

                # Map reads using bowtie
                
                echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${BASE} 21U RNAs with Bowtie..."
                
                bowtie -p 4 -r -S -v 3 -a --best --strata -m 400 /proj/ahmedlab/steve/seq/WS251/genome/bowtie/genome ./filtered/${BASE}_21u.txt ./bowtie_out/${BASE}_21u.sam
        
                # Convert to bam then sort

                echo $(date +"%m-%d-%Y_%H:%M")" Converting and sorting ${BASE}..."

                samtools view -bS ./bowtie_out/${BASE}_21u.sam | samtools sort -o ./bam/${BASE}_21u.bam -

                # Need to index the sorted bam files for visualization

                echo $(date +"%m-%d-%Y_%H:%M")" Indexing ${BASE}..."

                samtools index ./bam/${BASE}_21u.bam
        fi

        if [ "$FILTER_G" = false ] && [ "$FILTER_U" = false ]; then

                # Convert reads to raw format

                echo $(date +"%m-%d-%Y_%H:%M")" Extracting raw sequences from $file ..."
                
                python /proj/ahmedlab/steve/seq/util/small_rna_filter.py -o ./all_reads/reads/${BASE}.txt $file

                # Map reads using bowtie
                
                echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${BASE} with Bowtie..."
                
                bowtie -p 4 -r -S -v 3 -a --best --strata -m 400 /proj/ahmedlab/steve/seq/WS251/genome/bowtie/genome ./all_reads/reads/${BASE}.txt ./all_reads/sam/${BASE}.sam
        
                # Convert to bam then sort

                echo $(date +"%m-%d-%Y_%H:%M")" Converting and sorting ${BASE}..."

                samtools view -bS ./22g/sam/${BASE}.sam | samtools sort -o ./all_reads/bam/${BASE}.bam -

                # Need to index the sorted bam files for visualization

                echo $(date +"%m-%d-%Y_%H:%M")" Indexing ${BASE}..."

                samtools index ./all_reads/bam/${BASE}.bam
        fi

done

# Create count table using merge_counts.py

if [[ $COUNT == "true" ]]; then
        echo $(date +"%m-%d-%Y_%H:%M")"Merging count files into count table"
        python /proj/ahmedlab/steve/seq/util/merge_counts.py ./count
fi
