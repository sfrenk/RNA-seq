#!/usr/bin/env bash

# Hard variables

# Directory containing Snakemake and cluster.json files
snakedir='/nas/longleaf/home/sfrenk/pipelines/snakemake'

usage="Create directory with Snakemake files required for pipeline \n\n setup_dir -p <pipeline> -d <directory> \n\n pipelines: srna_telo, bowtie_srna, hisat2_stringtie"

pipeline=""

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

if [[ $pipeline == "" ]]; then
	echo "ERROR: invalid pipeline"
	exit 1
fi

# Determine pipeline file
case $pipeline in
	"srna_telo")
	snakefile="srna_telo"
	modules="anaconda python bowtie/1.1.2 samtools subread"
	;;
	"bowtie_srna")
	snakefile="bowtie_srna.Snakefile"
	modules="anaconda python bowtie/1.1.2 samtools subread"
	;;
	"hisat2_stringtie")
	snakefile='hisat2_stringtie.Snakefile'
	modules="python bbmap hisat2 samtools subread"
	;;
esac

# Copy over the snakefile
cp ${snakedir}/${snakefile} ./${snakefile}

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

# Create Snakmake command script

echo "module add $modules" > "run_snakemake.sh"
echo "snakemake -s $snakefile --cluster-config ${snakedir}/cluster.json -j 100 --cluster \"sbatch -n {cluster.n} -N {cluster.N} -t {cluster.time}\"" >> "run_snakemake.sh"
