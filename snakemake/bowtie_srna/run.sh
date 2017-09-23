#!/usr/bin/bash
#SBATCH -t 1-0
#SBATCH -N 1
#SBATCH -n 4

module add python anaconda bowtie samtools

snakemake --cores 4
