# Pipelines for performing RNA-seq analysis in C elegans

These pipelines have been written for use with the SLURM scheduling system, but they can be used on any system with a few slight modifications.

## Pipelines

### bowtie_sRNA.sh

Analysis pipeline for small RNA reads. Reads in raw/fasta/fastq format are mapped to specified reference sequence and can optionally be: 
	a) filtered for 22G and/or 21U RNAs
	b) counted
Note that read trimming and quality filtering are not included in this pipeline, and should therefore be performed using separate tools beforehand.

### hisat2.sh

Basic pipeline for mapping and counting single/paired end reads using hisat2, a much faster equivalent to tophat made by the same developers. This pipeline includes (optional) adapter trimming. Read counting can be done with subread and a reference transcript GTF file. Alternatively, stringtie can be used for de-novo transcriptome assembly.

## Script types

### Snakemake

I would recommend using the Snakemake scripts for running these pipelines. To run a pipeline, use setup_dir.sh to copy in the appropriate snakefile and create a bash script for running the pipeline.

setup_dir.sh -d (directory containing samples) -p (pipeline name)

### Bash scripts

Use these as a simple alternative to Snakemake (warning - these pipelines may not have been updated recently).
