#!/usr/bin/env bash

###############################################################################
# Basic pipeline for mapping and counting single/paired end reads using hisat2
###############################################################################

###############################################################################
# BEFORE RUNNING SCRIPT DO THE FOLLOWING:
###############################################################################

# Make sure modules are loaded:
#       trim_galore
#       hisat2
#       samtools
#       subread

# Check the location of the hard variables

###############################################################################
# Hard variables
###############################################################################

#       bowtie2 indices
index="/nas02/home/s/f/sfrenk/proj/seq/WS251/genome/bowtie2/genome"

#       gtf annotation file
gtf="/nas02/home/s/f/sfrenk/proj/seq/WS251/genes.gtf"

#       command to display software versions used during the run
modules=$(/nas02/apps/Modules/bin/modulecmd tcsh list 2>&1)

###############################################################################
###############################################################################


usage="
    USAGE
       step1:   load the following modules: trim_galore hisat2 samtools subread
       step2:   bash hisat2_genome.sh [options]  

    ARGUMENTS
        -d/--dir
        Directory containing read files (fastq.gz format).

        -p/--paired
        Use this option if fastq files contain paired-end reads. NOTE: if paired, each pair must consist of two files with the basename ending in '_r1' or '_r2' depending on respective orientation.

        -m/--multihits
        Maximum number of multiple hits allowed during hisat2 mapping (default = 1).
    "

# Set default parameters
paired=false
multihits=1


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
        -p|--paired)
        paired=true
        ;;
        --m|--multihits)
        multihits=$1
        shift
        ;;
    esac
shift
done

# Remove trailing "/" from input directory if present

if [[ ${dir:(-1)} == "/" ]]; then
    dir=${dir::${#dir}-1}
fi

# Print out loaded modules to keep a record of which software versions were used in this run

echo "$modules"

# module test

req_modules=("trim_galore" "hisat2" "samtools" "subread")

for i in ${req_modules[@]}; do
    if [[ $modules != *${i}* ]]; then
        echo "ERROR: Please load ${i}"
        exit 1
    fi
done

###############################################################################
###############################################################################

# Prepare directories

if [ ! -d "trimmed" ]; then
    mkdir trimmed
fi

if [ ! -d "hisat2_out" ]; then
    mkdir hisat2_out
fi

if [ ! -d "bam" ]; then
    mkdir bam
fi

if [ ! -d "count" ]; then
    mkdir count
fi

if [ -e "total_mapped_reads.txt" ]; then
    rm "total_mapped_reads.txt"
fi

echo "$(date +"%m-%d-%Y_%H:%M") Starting pipeline"

for file in ${dir}/*.fastq.gz; do
    
    skipfile=false

    if [[ $paired = true ]]; then
            
        # paired end

        if [[ ${file:(-11)} == "_1.fastq.gz" ]]; then
        
            FBASE=$(basename $file .fastq.gz)
            BASE=${FBASE%_1}

            echo $(date +"%m-%d-%Y_%H:%M")" Trimming ${BASE} with trim_galore..."

            trim_galore --dont_gzip -o ./trimmed --paired ${dir}/${BASE}_1.fastq.gz ${dir}/${BASE}_2.fastq.gz

            # Map reads using hisat2

            echo "$(date +"%m-%d-%Y_%H:%M") Mapping ${BASE} with hisat2... "        
            hisat2 --max-intronlen 12000 --no-mixed -k $multihits -p 4 -x ${index} -1 ./trimmed/${BASE}_1_val_1.fq -2 ./trimmed/${BASE}_2_val_2.fq -S ./hisat2_out/${BASE}.sam

        else

            # Avoid double mapping by skipping the r2 read file
                
            skipfile=true
        fi
    else

        # Single end

        BASE=$(basename $file .fastq.gz)

        # Trim reads

        echo $(date +"%m-%d-%Y_%H:%M")" Trimming ${BASE} with trim_galore..."

        trim_galore --dont_gzip -o ./trimmed ${dir}/${BASE}.fastq.gz

        # Map reads using hisat2

        echo "$(date +"%m-%d-%Y_%H:%M") Mapping ${BASE} with hisat2... "        
        hisat2 --max-intronlen 12000 --no-mixed -k $multihits -p 4 -x ${index} -U ./trimmed/${BASE}_trimmed.fq -S ./hisat2_out/${BASE}.sam
    fi

    if [[ $skipfile = false ]]; then

        echo $(date +"%m-%d-%Y_%H:%M")" Mapped ${BASE}"

        echo "$(date +"%m-%d-%Y_%H:%M") Sorting and indexing ${BASE}.bam"

        # Get rid of unmapped reads

        samtools view -h -F 4 ./hisat2_out/${BASE}.sam > ./bam/${BASE}.bam

        # Sort and index

        samtools sort -o ./bam/${BASE}_sorted.bam ./bam/${BASE}.bam

        samtools index ./bam/${BASE}_sorted.bam

        rm ./hisat2_out/${BASE}.sam
        rm ./bam/${BASE}.bam

        # Extract number of mapped reads

        total_mapped="$(samtools view -c ./bam/${base}_sorted.bam)"
        printf ${base}"\t"${total_mapped}"\n" >> total_mapped_reads.txt
    fi
done

echo $(date +"%m-%d-%Y_%H:%M")" Counting reads with featureCounts... "

# Count all files together so the counts will appear in one file

ARRAY=()

for file in ./bam/*_sorted.bam
do

ARRAY+=" "${file}

done

featureCounts -a ${gtf} -o ./count/counts.txt -T 4 -t exon -g transcript_id${ARRAY}
