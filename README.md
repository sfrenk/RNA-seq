# RNA-seq

## Pipelines for performing RNA-seq analysis in C elegans

These pipelines have been designed to automate RNA-seq data processing in a batch of samples. The final output includes sorted and indexed bam files and read count tables that can be further analyzed using the appropriate scripts provided in the seq-utilities repo. The pipelines also create a file named run_parameters.txt, which contains all the parameters and software versions used in the analysis.

The pipelines have been written for use on the Longleaf computing cluster at UNC - Chapel Hill so they are set up to be used with the SLURM scheduling system, but they can be used on any system with a few slight modifications.

## bowtie_sRNA.sh

Analysis pipeline for small RNA reads. Reads in raw/fasta/fastq format are mapped to specified reference sequence and can optionally be: 
	a) filtered for 22G and/or 21U RNAs
	b) counted
Note that read trimming and quality filtering are not included in this pipeline, and should therefore be performed using separate tools beforehand.

## hisat2.sh

Basic pipeline for mapping and counting single/paired end reads using hisat2, a much faster equivalent to tophat made by the same developers. This pipeline includes (optional) adapter trimming.
