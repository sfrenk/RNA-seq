#!/usr/bin/bash

# Hard variables

# Directory containing Snakemake and cluster.json files
snakefile_dir='/nas/longleaf/home/sfrenk/pipelines/snakemake'

usage="\nCreate directory with Snakemake files required for pipeline \n\n setup_dir -p <pipeline> -d <directory> \n\n pipelines: bowtie_srna, hisat2_rna, srna_telo\n\n"

pipeline=""

if [ -z "$1" ]; then
    printf "$usage"
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
	echo "ERROR: Please select pipeline: bowtie_srna or hisat2_rna"
	exit 1
fi

# Determine pipeline file
case $pipeline in
	"bowtie_srna"|"bowtie_sRNA")
	snakefile="bowtie_srna.Snakefile"
	;;
	"hisat2_rna"|"hisat2_RNA")
	snakefile='hisat2_rna.Snakefile'
	;;
	"srna_telo")
	snakefile="srna_telo.Snakefile"
	;;
	*)
	echo "ERROR: Invalid pipeline. Please select one of the following: bowtie_srna, hisat2_rna or srna_telo"
	exit 1
	;;
esac

# Copy over the snakefile
cp ${snakefile_dir}/${snakefile} ./${snakefile}

# Edit base directory in Snakefile
# Remove trailing "/" from dir if it's there 
input_dir="$(echo $dir |sed -r 's/\/$//')"
input_dir=\"${input_dir}\"
sed -i -e "s|^BASEDIR.*|BASEDIR = ${input_dir}|" $snakefile

# Determine file extension
extension="$(ls $dir | grep -Eo "\.[^/]+(\.gz)?$" | sort | uniq)"

# Check if there are multiple file extensions in the same directory
ext_count="$(ls $dir | grep -Eo "\.[^/]+(\.gz)?$" | sort | uniq | wc -l)"
if [[ $ext_count == 0 ]]; then
	echo "ERROR: Directory is empty!"
elif [[ $ext_count != 1 ]]; then
	echo "WARNING: Multiple file extensions found: using .fastq.gz"
	extension=".fastq.gz"
fi

# Edit extension and utils_dir in Snakefile
extension="\"${extension}\""
sed -i -e "s|^EXTENSION.*|EXTENSION = ${extension}|g" $snakefile
utils_dir="${snakefile_dir%/snakemake}/utils"
utils_dir="\"${utils_dir}\""
sed -i -e "s|^UTILS_DIR.*|UTILS_DIR = ${utils_dir}|g" $snakefile

# Create Snakmake command script
printf "#!/usr/bin/bash\n" > "run_snakemake.sh"
printf "#SBATCH -t 2-0\n\n" >> "run_snakemake.sh"
printf "module add python\n\n" >> "run_snakemake.sh"
printf "snakemake -s $snakefile --keep-going --rerun-incomplete --cluster-config ${snakefile_dir}/cluster.json -j 100 --cluster \"sbatch -n {cluster.n} -N {cluster.N} -t {cluster.time}\"\n" >> run_snakemake.sh
