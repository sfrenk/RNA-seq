# Pipelines for performing RNA-seq analysis in C elegans

These pipelines have been written for use with the SLURM scheduling system, but they can be used on any system with a few slight modifications.

## Installation

```
git clone https://github.com/sfrenk/sTeloMap.
cd rna-seq_pipelines/
bash install.sh
```


## Usage

Once the repo has been installed, running the pipeline is simple:

1. Use [setup_dir.sh](utils/setup_dir.sh) (located in the utils subdirectory) to set up the working directory for analysis:

```
path/to/utils/setup_dir.sh -p <pipeline> -d <directory containing reads>
```

2. Run the snakemake script

```
bash run_snakemake.sh
```

## Pipelines

### [bowtie_sRNA](snakemake/bowtie_srna.Snakefile)

Analysis pipeline for small RNA reads. Reads in raw/fasta/fastq format are mapped to specified reference sequence and can optionally be: 
	a) filtered for 22G and/or 21U RNAs
	b) counted
Note that read trimming and quality filtering are not included in this pipeline, and should therefore be performed using separate tools beforehand.

### [hisat2_RNA](snakemake/hisat2_rna.Snakefile)

Basic pipeline for mapping and counting single/paired end reads using hisat2, a much faster equivalent to tophat made by the same developers. This pipeline includes (optional) adapter trimming. Read counting can be done with subread and a reference transcript GTF file. Alternatively, stringtie can be used for de-novo transcriptome assembly.


## Downstream analysis

The final output file for hisat2_RNA (and bowtie_sRNA if using the genome reference) is a table produced by subread called "counts.txt", located in the count subdirectory of the working directory. The utils directory contains the script [DESeq2.R](utils/DESeq2.R) which can be used to perform differential gene expression analysis using the count data:

```
Rscript /path/to/utils/DESeq2.R -o <output filename> -c <"control" samples> -t <"treatment" samples> <counts.txt>
```

"control" and "treatment" are just names to distinguish the two groups of samples that are being compared. Note that the log2-fold change values represent the change (increase or decrease) in expression in the treatment compared with the control group.
