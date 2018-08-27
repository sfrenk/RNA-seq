#!/usr/bin/env bash

# Hard variables

# Directory containing Snakemake and cluster.json files
snakedir='/nas/longleaf/home/sfrenk/pipelines/snakemake'

usage="Create directory with Snakemake files required for pipeline \n\n setup_dir -p <pipeline> -d <directory> \n\n pipelines: srna_telo, bowtie_srna, hisat2_stringtie, chip_seq, gatk"

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
		printf "$usage"
		exit
		;;
    esac
    shift
done


if [[ ! -d $dir ]]; then
	echo "ERROR: Invalid directory"
	exit 1
fi

if [[ $pipeline == "" ]]; then
	echo "ERROR: Please select pipeline"
	exit 1
fi

# Determine pipeline file
case $pipeline in
	"srna_telo")
	snakefile="srna_telo.Snakefile"
	;;
	"bowtie_srna")
	snakefile="bowtie_srna.Snakefile"
	;;
	"hisat2_stringtie")
	snakefile='hisat2_stringtie.Snakefile'
	;;
	"chip_seq")
	snakefile='chip_seq.Snakefile'
	;;
	"gatk"|"call_variants")
	snakefile='call_variants.Snakefile'
	;;
	"gatk_human"|"call_variants_human")
	snakefile='call_variants_human_germline.Snakefile'
	modules="python"
	;;
	*)
	echo "ERROR: Invalid pipeline. Please select one of the following: bowtie_srna, hisat2_stringtie, chip_seq, gatk"
	exit 1
	;;
esac

# Copy over the snakefile
cp ${snakedir}/${snakefile} ./${snakefile}

# Edit base directory in Snakefile
# Remove trailing "/" from dir if it's there 
dir_name="$(echo $dir |sed -r 's/\/$//')"
sed -r -i -e "s,^BASEDIR.*,BASEDIR = \"${dir_name}\"," "$snakefile"

# Determine file extension
extension="$(ls $dir | grep -Eo "\.[^/]+(\.gz)?$" | sort | uniq)"

# Check if there are multiple file extensions in the same directory
ext_count="$(echo $extension | wc -l)"

if [[ $ext_count == 0 ]]; then
	echo "ERROR: Directory is empty!"
elif [[ $ext_count -gt 1 ]]; then
	echo "WARNING: Multiple file extensions found: using .fastq.gz"
	extension=".fastq.gz"
fi

# Edit extension in Snakefile
extension="\"${extension}\""
sed -i -r -e "s/^EXTENSION.*/EXTENSION = ${extension}/g" "$snakefile"

# Create Snakmake command script
printf "#!/usr/bin/bash\n" > "run_snakemake.sh"
printf "#SBATCH -t 2-0\n\n" >> "run_snakemake.sh"
printf "module add python\n\n" >> "run_snakemake.sh"
printf "snakemake -s $snakefile --keep-going --rerun-incomplete --cluster-config ${snakedir}/cluster.json -j 100 --cluster \"sbatch -n {cluster.n} -N {cluster.N} -t {cluster.time}\"\n" >> "run_snakemake.sh"
