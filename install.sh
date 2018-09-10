#!/usr/bin/bash -e

echo "Setting up RNA-seq pipelines"

# Find snakefile and utils directories
script_dir="$(pwd)"
snakefile_dir="${script_dir}/snakemake"
utils_dir="${snakefile_dir%/snakemake}/utils"
utils_dir="\"${utils_dir}\""

# Edit setup_dir.sh to contain correct snakefile directory
sed -i -e "s|^snakefile_dir.*|snakefile_dir=${snakefile_dir}|g" "${snakefile_dir}/setup_dir.sh"

# Edit Snakefiles to have the correct utils directory
snakefiles=("${snakefile_dir}/*.Snakefile")

for i in snakefiles; do
	sed -i -e "s|^UTILS_DIR.*|UTILS_DIR = ${utils_dir}|g" $i
done

echo "Done!"
