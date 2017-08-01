#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 1-0

module add bowtie samtools python/2.7.12 subread
module list

###############################################################################
# bowtie_sRNA.sh: Analysis pipeline for small RNA reads (version 2.0)
###############################################################################

# This script maps reads in raw, fasta or fastq format to a specified reference sequence and can also perform the following tasks:
#       filtering for 22G and/or 21U RNAs
#       counting the number of reads mapping to each feature

# NOTE: Quality filtering / adaptor trimming is NOT covered in the pipeline

# The location of the following files may have to be modified in this script:
# bowtie index files
# merge_counts.py and small_rna_filter.py

###############################################################################
###############################################################################

usage="
    This script maps reads in raw, fasta or fastq format to a specified reference sequence and can also perform the following tasks:
        filtering for 22G and/or 21U RNAs
        counting the number of reads mapping to each feature

    This script should be run with four processors

    USAGE
       
       bash bowtie_sRNA.sh [options]  

    ARGUMENTS
        -d/--dir
        directory containing read files in fasta.gz, txt.gz or fastq.gz format
        
        -r/--reference
        Mapping reference: genome, transposons, mrna, mirna, pirna, telomere, rdna or e_coli (default = genome)

        -n/--no_filter
        Do not filter reads. Use this option for non-small RNA libraries. Files must be in .fastq.gz format.

        -f/--filter 
        Filter reads based on the first nucleotide (eg. to select 22g RNAs, use options -f g -s 22,22)

        -s/--size
        Specify size range of reads to keep by providing min and max length seperated by a comma (eg. to keep reads between 19 and 24 nucleotides (inclusive), use '-s 19,24' (Default size range: 18 to 30 nucleotides)

        -t/--trim_a
        Trim 3' A nucleotides from reads before filtering/mapping

        -c/--count
        count reads and compile read counts from each sample into a count table

        -m/--mismatch
        set maximum number of base mismatches allowed during mapping (default = 0)

        -a/--antisense
        if --count option is selected, only reads that map antisense to reference sequence are kept

        -l/--all
        If there are more than one best alignments for a read, report all alignments (by default the read is assigned to one location at random)
    "
# Set default parameters

ref="genome"
mismatch=0
filter="A,T,C,G"
size="18,30"
count=false
antisense=false
trim_a=""
no_filter=false
all=false

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
            -m|--mismatch)
            mismatch="$2"
            shift
            ;;
            -f|--filter)
            filter="$2"
            shift
            ;;
            -s|--size)
            size="$2"
            shift
            ;;
            -n|--no_filter)
            no_filter=true
            ;;
            -c|--count)
            count=true
            ;;
            -a|--antisense)
            antisense=true
            ;;
            -t|--trim_a)
            trim_a="-a"
            ;;
            -l|--all)
            all=true
            ;;
    esac
    shift
done

# Remove trailing "/" from dir directory if present

if [[ ${dir:(-1)} == "/" ]]; then
    dir=${dir::${#dir}-1}
fi

# Select bowtie index based on reference

case $ref in 
    "mrna")
    index="/nas02/home/s/f/sfrenk/proj/seq/WS251/mrna/bowtie/mrna"
    ;;
    "transposons")
    index="/proj/ahmedlab/steve/seq/transposons/bowtie/transposon"
    ;;
    "genome")
    index="/proj/ahmedlab/steve/seq/WS251/genome/bowtie/genome"
    ;;
    "mirna")
    index="/proj/ahmedlab/steve/seq/mirna/bowtie/mirna"
    ;;
    "rdna")
    index="/proj/ahmedlab/steve/seq/rdna/bowtie/rdna"
    ;;
    "telomere")
    index="/proj/ahmedlab/steve/seq/telomere/bowtie/telomere"
    ;;
    "trna")
    index="/nas02/home/s/f/sfrenk/seq/trna/bowtie/trna"
    ;;
    "pirna")
    index="/nas02/home/s/f/sfrenk/seq/pirna/bowtie/pirna"
    ;;
    "e_coli")
    index="/proj/ahmedlab/steve/seq/e_coli_b/genome/bowtie/genome"
esac

# Print run parameters to file

if [[ ! -f run_parameters.txt ]]; then
    printf $(date +"%m-%d-%Y_%H:%M")"\n\nPipeline: bowtie small RNA\n\nParameters:\n\tsample directory: ${dir}\n\tref: ${ref}\n\tmismatches: ${mismatch}\n\tfirst nucleotide: ${filter}\n\tsize range: ${size}\n\tcount reads: ${count}\n\tcount only antisense reads: ${antisense}\n\treport all alignments: ${all}\n\nM\n\nSamples:" > run_parameters.txt

    module list &>> run_parameters.txt
    printf "\n" >> run_parameters.txt
fi

###############################################################################
###############################################################################

# Make directories

if [ ! -d "filtered" ]; then
    mkdir filtered
fi
if [ ! -d "bowtie_out" ]; then
    mkdir bowtie_out
fi
if [ ! -d "count" ] && [[ ${count} = true ]]; then
    mkdir count
fi
if [ ! -d "bam" ]; then
    mkdir bam
fi


# Start pipeline

echo $(date +"%m-%d-%Y_%H:%M")" Starting pipeline..."

shopt -s nullglob

files=(${dir}/*.fastq.gz)

for file in ${files[@]}; do
        
    # Define basename of sample based on file extension
    
    if [[ ${file:(-9)} == ".fastq.gz" ]]; then
        base=$(basename $file .fastq.gz)
    elif [[ ${file:(-6)} == ".fq.gz" ]]; then
        base=$(basename $file .fq.gz)
    elif [[ ${file:(-7)} == ".txt.gz" ]]; then
        base=$(basename $file .txt.gz)
    elif [[ ${file:(-9)} == ".fasta.gz" ]]; then
        base=$(basename $file .fasta.gz)
    elif [[ ${file:(-6)} == ".fa.gz" ]]; then
        base=$(basename $file .fa.gz)
    fi

    echo $(date +"%m-%d-%Y_%H:%M")" ################ processing ${base} ################"

    printf "\n\t"$base >> run_parameters.txt


    # Extract 22G and or 21U RNAs or convert all reads to raw format
    
    python /proj/ahmedlab/steve/seq/util/small_rna_filter.py -f $filter -s $size -o ./filtered/${base}.txt $trim_a $file

    # Map reads using Bowtie
    
    echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${base} with Bowtie..."

    if [[ $all = true ]]; then
        
        bowtie --best --strata -a -r -S -v $mismatch -p 8 $index ./filtered/${base}.txt ./bowtie_out/${base}.sam
    else
        
        bowtie --best --strata -m 1 -r -S -v $mismatch -p 8 $index ./filtered/${base}.txt ./bowtie_out/${base}.sam
    fi


    echo $(date +"%m-%d-%Y_%H:%M")" Mapped ${base}"
        
    # Convert to bam then sort

    echo $(date +"%m-%d-%Y_%H:%M")" Converting and sorting ${base}..."

    samtools view -bh -F 4 ./bowtie_out/${base}.sam > ./bam/${base}.bam

    samtools sort -o ./bam/${base}_sorted.bam ./bam/${base}.bam

    rm ./bam/${base}.bam

    # Need to index the sorted bam files for visualization

    echo $(date +"%m-%d-%Y_%H:%M")" indexing ${base}..."

    samtools index ./bam/${base}_sorted.bam


    # Count reads by location/feature
    
    if [[ $count = true ]] && [[ $ref != "genome" ]] && [[ $ref != "e_coli" ]]; then
        
        if [[ $antisense = false ]]; then
            echo $(date +"%m-%d-%Y_%H:%M")" counting reads"
            
            # Extract and count reads
            # The 260 flag marks unmapped reads/secondary mappings
            
            awk -F$'\t' '$2 != "260" ' ./bowtie_out/${base}.sam | grep -v @ | cut -f 3 | sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' > ./count/${base}_counts.txt
            
        else
            
            # Only count reads that map antisense to features
            
            echo $(date +"%m-%d-%Y_%H:%M")" counting antisense reads"

            awk -F$'\t' '$2 == "16" ' ./bowtie_out/${base}.sam | grep -v @ | cut -f 3 |  sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' > ./count/${base}_counts.txt
        fi
    fi

    rm ./bowtie_out/${base}.sam

    # Get total number of mapped reads

    total_mapped="$(samtools view -c ./bam/${base}_sorted.bam)"
    printf ${base}"\t"${total_mapped}"\n" >> total_mapped_reads.txt

    echo $(date +"%m-%d-%Y_%H:%M")" ################ ${base} has been processed ################"
done

# Create count table using merge_counts.py

if [[ $count = true ]]; then

    if [[ $ref == "genome" ]] || [[ $ref == "e_coli" ]]; then

        # Set up specific parameters for C elegans or E coli
        if [[ $ref == "genome" ]]; then
        
            gtf="/nas02/home/s/f/sfrenk/proj/seq/WS251/genes.gtf"
            feature="exon"
            att="gene_name"

        else

            gtf="/proj/ahmedlab/steve/seq/e_coli_b/annotation/genes.gff"
            feature="gene"
            att="ID"

        fi

        ARRAY=()

        for file in ./bam/*_sorted.bam
        do

        ARRAY+=" "${file}

        done

        if [[ $antisense = true ]]; then
            
            featureCounts -M -s 2 -a ${gtf} -o ./count/counts_output.txt -T 4 -t $feature -g ${att}${ARRAY}

        else
            featureCounts -M -a ${gtf} -o ./count/counts_output.txt -T 4 -t $feature -g ${att}${ARRAY}
        fi
        
        # Clean up sample names in count file

        sed -e "2s/\.\/bam\///g" ./count/counts_output.txt | sed -e "2s/_sorted\.bam//g" > ./count/counts.txt

    else

        echo $(date +"%m-%d-%Y_%H:%M")" Merging count files into count table"
        python /proj/ahmedlab/steve/seq/util/merge_counts.py ./count

    fi
fi
