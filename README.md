# Pipelines for performing RNA-seq analysis in C elegans

These pipelines have been written for use with the SLURM scheduling system, but they can be used on any system with a few slight modifications.

## Installation

```
git clone https://github.com/sfrenk/sTeloMap.git

bash install.sh
```


## Usage

Once the repo has been installed, running the pipeline is simple:

1. Use setup_working_dir.sh to set up the working directory for analysis:

```
setup_dir.sh -p <pipeline> -d <directory containing reads>
```

2. Run the snakemake script

```
bash run_snakemake.sh
```


## Pipelines

### bowtie_sRNA

Analysis pipeline for small RNA reads. Reads in raw/fasta/fastq format are mapped to specified reference sequence and can optionally be: 
	a) filtered for 22G and/or 21U RNAs
	b) counted
Note that read trimming and quality filtering are not included in this pipeline, and should therefore be performed using separate tools beforehand.

### hisat2_RNA

Basic pipeline for mapping and counting single/paired end reads using hisat2, a much faster equivalent to tophat made by the same developers. This pipeline includes (optional) adapter trimming. Read counting can be done with subread and a reference transcript GTF file. Alternatively, stringtie can be used for de-novo transcriptome assembly.
