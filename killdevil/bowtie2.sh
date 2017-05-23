#!/usr/bin/env bash

###############################################################################
# Basic pipeline for mapping and counting single/paired end reads using bowtie2
###############################################################################

###############################################################################
# Hard variables
###############################################################################

#       gtf annotation file
gtf="/nas02/home/s/f/sfrenk/proj/seq/WS251/genes.gtf"

#       command to display software versions used during the run
modules=$(/nas02/apps/Modules/bin/modulecmd tcsh list 2>&1)

###############################################################################
###############################################################################


usage="
    Basic bowtie2 pipeline for RNA-seq analysis. Note that the --very-sensitive option is used during bowtie2 alignment.

    This script should be run with four processors

    USAGE
       step1:   load the following modules: trim_galore bowtie2 samtools subread
       step2:   bash bowtie2_genome.sh [options]  

    ARGUMENTS
        -d/--dir
        Directory containing read files (fastq.gz format).

        -p/--paired
        Use this option if fastq files contain paired-end reads. NOTE: if paired, each pair must consist of two files with the basename ending in '_r1' or '_r2' depending on respective orientation.

        -a/--all
        Report all alignments (by default, this script calls the bowtie2 default reporting option, which is to report the best alignment at random)

        -t/--trim
        Trim/filter reads with trim_galore before mapping

        -l/--local
        Perform local (rather than end-to-end) alignment

        -r/--reference
        Mapping reference: genome, transposons, mrna, mirna, pirna, telomere or rdna (default = genome)

        -a/--antisense
        if --count option is selected, only reads that map antisense to reference sequence are kept
    "

# Set default parameters
paired=false
all=""
trim=false
alignment_mode="--very-sensitive"
ref="genome"
antisense=false

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
        --a|--all)
        all=" -a"
        ;;
        -t|--trim)
        trim=true
        ;;
        -l|--local)
        alignment_mode="--very-sensitive-local"
        ;;
        -r|--ref)
        ref="$2"
        shift
        ;;
        -a|--antisense)
        antisense=true
        ;;
    esac
shift
done

# Remove trailing "/" from input directory if present

if [[ ${dir:(-1)} == "/" ]]; then
    dir=${dir::${#dir}-1}
fi

# Select bowtie index based on reference

case $ref in 
    "mrna")
    index="/nas02/home/s/f/sfrenk/proj/seq/WS251/mrna/bowtie2/mrna"
    ;;
    "transposons")
    index="/proj/ahmedlab/steve/seq/transposons/bowtie2/transposon"
    ;;
    "genome")
    index="/proj/ahmedlab/steve/seq/WS251/genome/bowtie2/genome"
    ;;
    "mirna")
    index="/proj/ahmedlab/steve/seq/mirna/bowtie2/mirna"
    ;;
    "rdna")
    index="/proj/ahmedlab/steve/seq/rdna/bowtie2/rdna"
    ;;
    "telomere")
    index="/proj/ahmedlab/steve/seq/telomere/bowtie2/telomere"
    ;;
    "trna")
    index="/nas02/home/s/f/sfrenk/seq/trna/bowtie2/trna"
    ;;
    "pirna")
    index="/nas02/home/s/f/sfrenk/seq/pirna/bowtie2/pirna"
    ;;
esac

# Print out loaded modules to keep a record of which software versions were used in this run

echo "$modules"

# module test

req_modules=("bowtie2" "samtools")

for i in ${req_modules[@]}; do
    if [[ "$modules" != *${i}* ]]; then
        echo "ERROR: Please load ${i}"
        exit 1
    fi
done

if [[ $trim = true ]] && [[ "$modules" != *"trim_galore"* ]]; then
    echo "ERROR: Please load trim_galore"
    exit 1
fi

if [[ $ref == "genome" ]] && [[ "$modules" != *"subread"* ]]; then
    echo "ERROR: Please load subread"
    exit 1
fi

if [[ $ref != "genome" ]] && [[ "$modules" != *"python"* ]]; then
    echo "ERROR: Please load python"
    exit 1
fi


# Print run parameters to file

if [ -e "run_parameters.txt" ]; then
    rm "run_parameters.txt"
fi

printf $(date +"%m-%d-%Y_%H:%M")"\n\nPipeline: bowtie2\n\nParameters:sample directory: ${dir}\n\treport all alignments?: ${all}\n\tpaired end: ${paired}\n\ttrim: ${trim}\n\nalignment mode: ${alignment_mode}\n\nreference: ${ref}\n\ncount antisense only?: ${antisense}\n\nModules: ${modules}\n\nSamples:" > run_parameters.txt

###############################################################################
###############################################################################

# Prepare directories

if [ ! -d "trimmed" ] && [[ $trim = true ]]; then
    mkdir trimmed
fi

if [ ! -d "bowtie2_out" ]; then
    mkdir bowtie2_out
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
        
            Fbase=$(basename $file .fastq.gz)
            base=${Fbase%_1}

            printf "\n\t"$base >> run_parameters.txt

            echo $(date +"%m-%d-%Y_%H:%M")" Trimming ${base} with trim_galore..."

            if [[ $trim = true ]]; then

                # Trim reads

                echo $(date +"%m-%d-%Y_%H:%M")" Trimming ${base} with trim_galore..."
                trim_galore --dont_gzip -o ./trimmed --paired ${dir}/${base}_1.fastq.gz ${dir}/${base}_2.fastq.gz

                fastq_r1="./trimmed/${base}_1_val_1.fq"
                fastq_r2="/trimmed/${base}_2_val_2.fq"
            else

                fastq_r1="${dir}/${base}_1.fastq.gz"
                fastq_r2="${dir}/${base}_2.fastq.gz"
            fi

            # Map reads using bowtie2

            echo "$(date +"%m-%d-%Y_%H:%M") Mapping ${base} with bowtie2... "        
            bowtie2 -p 4 $alignment_mode --no-mixed${all} -x ${index} -1 $fastq_r1 -2 $fastq_r2 -S ./bowtie2_out/${base}.sam

        else

            # Avoid double mapping by skipping the r2 read file
                
            skipfile=true
        fi
    else

        # Single end

        base=$(basename $file .fastq.gz)

        printf "\n\t"$base >> run_parameters.txt

        if [[ $trim = true ]]; then

            # Trim reads

            echo $(date +"%m-%d-%Y_%H:%M")" Trimming ${base} with trim_galore..."

            trim_galore --dont_gzip -o ./trimmed ${dir}/${base}.fastq.gz
            fastq_file="./trimmed/${base}_trimmed.fq"
        else

            fastq_file="${dir}/${base}.fastq.gz"
        fi

        # Map reads using bowtie2

        echo "$(date +"%m-%d-%Y_%H:%M") Mapping ${base} with bowtie2... "        
        bowtie2 -p 4 ${alignment_mode}${all} -x ${index} -U $fastq_file -S ./bowtie2_out/${base}.sam
    fi

    if [[ $skipfile = false ]]; then

        echo $(date +"%m-%d-%Y_%H:%M")" Mapped ${base}"

        echo "$(date +"%m-%d-%Y_%H:%M") Sorting and indexing ${base}.bam"

        # Get rid of unmapped reads

        samtools view -h -F 4 ./bowtie2_out/${base}.sam > ./bam/${base}.bam

        # Sort and index

        samtools sort -o ./bam/${base}_sorted.bam ./bam/${base}.bam

        samtools index ./bam/${base}_sorted.bam

        rm ./bam/${base}.bam

        # Extract number of mapped reads

        total_mapped="$(samtools view -c ./bam/${base}_sorted.bam)"
        printf ${base}"\t"${total_mapped}"\n" >> total_mapped_reads.txt
    fi

    # Count reads by location/feature
    
    if [[ $ref != "genome" ]]; then
        
        if [[ $antisense = false ]]; then
            echo $(date +"%m-%d-%Y_%H:%M")" counting reads"
            
            # Extract and count reads
            # The 260 flag marks unmapped reads/secondary mappings
            
            awk -F$'\t' '$2 != "260" ' ./bowtie2_out/${base}.sam | grep -v @ | cut -f 3 | sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' > ./count/${base}_counts.txt
            
        else
            
            # Only count reads that map antisense to features
            
            echo $(date +"%m-%d-%Y_%H:%M")" counting antisense reads"

            awk -F$'\t' '$2 == "16" ' ./bowtie2_out/${base}.sam | grep -v @ | cut -f 3 |  sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' > ./count/${base}_counts.txt
        fi

        rm ./bowtie2_out/${base}.sam
    fi
done

echo $(date +"%m-%d-%Y_%H:%M")" Counting reads with featureCounts... "

# Count all files together so the counts will appear in one file
if [[ $ref == "genome" ]]; then
    ARRAY=()

    for file in ./bam/*_sorted.bam
    do

        ARRAY+=" "${file}

    done

    featureCounts -a ${gtf} -o ./count/counts.txt -T 4 -t exon -g gene_name${ARRAY}
fi
