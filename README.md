# RNA-seq

## Pipelines for performing RNA-seq analysis in C elegans

## bowtie2_pipeline.sh

Basic gene expression analysis pipeline for mapping single or paired end fastq reads to mRNA or transposon consensus sequences with bowtie2.

## bowtie_sRNA.sh

Analysis pipeline for small RNA reads. Reads in raw/fasta/fastq format are mapped to specified reference sequence and can optionally be: 
	a) filtered for 22G and/or 21U RNAs
	b) counted

## tophat_genome.sh

Basic pipeline for mapping and counting single/paired end reads using tophat. This pipeline has been superseded by the hisat2 pipelines.

## hisat2_genome.sh

Basic pipeline for mapping and counting single/paired end reads using hisat2, a much faster equivalent to tophat made by the same developers.

## hisat2_fasta.sh

A hisat2 pipeline that takes reads in fasta format and does not include a trimming/filtering step. This pipeline is largely used for mapping small RNA reads to the genome.
