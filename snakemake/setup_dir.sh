#!/usr/bin/env bash

usage=
"Create directory with Snakemake files required for pipeline

setup_dir -p <pipeline> -d <directory>

pipelines: srna_telo (default), bowtie_srna"

pipeline="srna_telo"

if [ -z "$1" ]; then
    echo "$usage"
    exit
fi

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
    	-p|--pipeline)
		pipeline="$2"
		shift
		;;
        -d|--dir)
        dir="$2"
        shift
        ;;
    esac
    shift
done


if [[ ! -d $dir ]]; then
	echo "ERROR: invalid directory"
	exit 1
fi

# Determine pipeline dir

case $pipeline in
	"srna_telo")
	script_dir="/nas/longleaf/home/sfrenk/scripts/charlie_paper/snakemake"
	;;
	"bowtie_srna")
	script_dir="/nas/longleaf/home/sfrenk/pipelines/snakemake/bowtie_srna"
	;;
	"hisat2_stringtie")
	script_dir='/nas/longleaf/home/sfrenk/pipelines/snakemake/hisat2_stringtie'
	;;
esac

# Copy over the necessary files
cp "${script_dir}/Snakefile" .
cp "${script_dir}/cluster.json" .
cp "${script_dir}/config.json" .

# Edit base directory in Config file
base="$(basename ${dir})"
sed -i -r -e "s,\"basedir\": \".*\",\"basedir\": \"${srna_dir}\"," "./${base}/config.json"

# Determine file extension
extension="$(ls $srna_dir | grep -Eo "\.[^/]+" | sort | uniq)"

if [[ "${#extension[@]}" != 1 ]]; then
	echo "ERROR: All files must have the same extension"

else

	# Edit extension in config file
	sed -i -r -e "s,\"extension\": \".*\",\"extension\": \"${extension}\"," "./${base}/config.json"
fi
