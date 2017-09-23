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
        Mapping reference: genome, transposons, mrna, mirna, pirna, telomere, rdna, csr1_targets or e_coli (default = genome)

        -q/--fastq_mode
        Output reads in fastq format after filtering. By default, reads are converted to raw format, which removes quality information 

        -f/--filter 
        Filter reads based on the first nucleotide (eg. to select 22g RNAs, use options -f g -s 22,22)

        -s/--size
        Specify size range of reads to keep by providing min and max length seperated by a comma (eg. to keep reads between 19 and 24 nucleotides (inclusive), use '-s 19,24' (Default size range: 18 to 30 nucleotides)

        -t/--trim
        Trim 3' nucleotides of this base from 3' ends of reads before filtering/mapping (eg. use -t A to trim any 3' A nucleotides)

        -c/--count
        Count reads and compile read counts from each sample into a count table

        -m/--mismatch
        set maximum number of base mismatches allowed during mapping (default = 0)

        -a/--antisense
        If --count option is selected, only reads that map antisense to reference sequence are kept

        -l/--multi
        How to deal with reads with more than one best alignment:

            random (default: randomly assign the read to one of the locations
            discard: discard the read
            all: report all alignments

        --min_trim_length
        Stop trimming 3' nucleotides when read has this many nucleotides

        --csr1
        Use the following procedure to detect polyudyrilated csr-1 sRNAs:
            1. Map reads
            2. Take all unmapped reads ending with T and remove the final T
            3. Map the T-trimmed reads
            4. Keep repeating 1-3 until there are no T-trimmed unmapped reads of >= min size left 
    "
# Set default parameters

ref="genome"
mismatch=0
filter="A,T,C,G"
size="18,30"
count=false
antisense=false
trim=""
fastq_mode=false
raw_flag="-r"
multi_option="random"
csr=false
min_trim_length=""

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
        -q|--fastq)
        fastq_mode=true
        ;;
        -c|--count)
        count=true
        ;;
        -a|--antisense)
        antisense=true
        ;;
        -t|--trim)
        trim="-t $2"
        ;;
        -l|--multi)
        multi_option="$2"
        shift
        ;;
        --min_trim_length)
        min_trim_length="-m $2"
        shift
        ;;
        --csr1)
        csr=true
        # Need to use fastq mode for csr1 mode so that unmapped reads file doesn't look weird
        fastq_mode=true
        ;;
    esac
    shift
done

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
    "csr1")
    index="/nas/longleaf/home/sfrenk/proj/seq/WS251/srna/csr1_targets/bowtie/csr1_targets"
    ;;
    "e_coli")
    index="/proj/ahmedlab/steve/seq/e_coli_b/genome/bowtie/genome"
    ;;
esac

# parse multi_option

valid_multi_options=("discard" "random" "all")

if ! [[ "${valid_multi_options[@]}" =~ "$multi_option" ]]; then
    echo "ERROR: Please supply valid -l/--multi_option argument"
    exit
fi

case $multi_option in
    "discard")
    multi_flag="-m 1"
    ;;
    "random")
    multi_flag="-M 1"
    ;;
    "all")
    multi_flag="-a"
    ;;
esac

# Print run parameters to file

if [[ ! -f run_parameters.txt ]]; then
    printf $(date +"%m-%d-%Y_%H:%M")"\n\nPipeline: bowtie small RNA\n\nParameters:\n\tsample directory: ${dir}\n\tref: ${ref}\n\tmismatches: ${mismatch}\n\tfirst nucleotide: ${filter}\n\tsize range: ${size}\n\tcount reads: ${count}\n\tcount only antisense reads: ${antisense}\n\tmulti_option: ${multi_option}\n\tcsr1 mode (T trim): ${csr}\n\n" > run_parameters.txt

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
if [ ! -d "count" ] && [[ $count = true ]]; then
    mkdir count
fi
if [ ! -d "bam" ]; then
    mkdir bam
fi
if [ ! -d "fastq" ] && [[ $csr = true ]]; then
    mkdir fastq
fi

# Start pipeline

echo $(date +"%m-%d-%Y_%H:%M")" Starting pipeline..."

shopt -s nullglob

files=(${dir}/*.gz)

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

    printf "${base}\n" >> run_parameters.txt


    if [[ $fastq_mode = false ]]; then
    
        # Extract 22G and or 21U RNAs or convert all reads to raw format

        small_rna_filter -f $filter -s $size -o ./filtered/${base}.txt $trim $min_trim_length $file

        reads_file=./filtered/${base}.txt
    
    else

        # Need to switch off the -r argument in bowtie

        raw_flag="-q"
        small_rna_filter -f $filter -s $size -o ./filtered/${base}.fastq $trim $min_trim_length $file

        reads_file=./filtered/${base}.fastq
    fi

    # Map reads using Bowtie
    
    if [[ $csr = true ]]; then
        t_trim=1
        counter=0

        while [[ $t_trim -gt 0 && $count -le 10 ]]; do

            echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${base} in CSR-1 mode..."

            bowtie --best --strata $raw_flag $multi_flag -S -v $mismatch -p $SLURM_NTASKS --un ./fastq/${base}_unmapped.fastq $index $reads_file ./bowtie_out/${base}_${counter}.sam
            
            if [[ -f ./fastq/${base}_unmapped.fastq ]]; then
                counter=$((counter+1))

                small_rna_filter -f $filter -s $size -u T -o ./fastq/${base}_t_${counter}.fastq ./fastq/${base}_unmapped.fastq

                reads_file=./fastq/${base}_t_${counter}.fastq
                t_trim=$(wc -l < $reads_file)
            
            else
                
                t_trim=0
            fi

            rm ./fastq/${base}_unmapped.fastq
        
        done

        # Merge sam files
        
        header=./bowtie_out/${base}_0.sam
        files=(./bowtie_out/${base}_*.sam)

        (grep ^@ $header; for f in $files; do grep -v ^@ $f; done) > ./bowtie_out/${base}.sam

        rm ./bowtie_out/${base}_*.sam

    else

        echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${base} with Bowtie..."
        
        bowtie --best --strata $raw_flag $multi_flag -S -v $mismatch -p $SLURM_NTASKS $index $reads_file ./bowtie_out/${base}.sam
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
