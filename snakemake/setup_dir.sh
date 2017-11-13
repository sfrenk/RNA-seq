#!/usr/bin/env bash

# Hard variables
snakedir='/nas/longleaf/home/sfrenk/pipelines/snakemake'

usage="Create directory with Snakemake files required for pipeline \n\n setup_dir -p <pipeline> -d <directory> \n\n pipelines: srna_telo (default), bowtie_srna"

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
        -h|--help)
		echo "$usage"
		exit
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
	snakefile=""
	;;
	"bowtie_srna")
	snakefile="bowtie_srna.Snakefile"
	;;
	"hisat2_stringtie")
	snakefile='hisat2_stringtie.Snakefile'
	;;
esac

# Copy over the necessary files
cp ${snakedir}/${snakefile} ./${snakefile}
cp ${snakedir}/cluster.json .

# Edit base directory in Snakefile
base="$(basename ${dir})"
sed -r -i -e "s,^BASEDIR.*,BASEDIR = \"${dir}\"," "$snakefile"

# Determine file extension
extension="$(ls $dir | grep -Eo "\.[^/]+" | sort | uniq)"

if [[ "${#extension[@]}" != 1 ]]; then
	echo "ERROR: All files must have the same extension"

else

	# Edit extension in Snakefile
	extension="\"${extension}\""
	sed -i -r -e "s/^EXTENSION.*/EXTENSION = ${extension}/g" "$snakefile"
fi
