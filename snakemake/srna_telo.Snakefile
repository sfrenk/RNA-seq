import glob
import re
import sys
import os

### NEW AND IMPROVED PIPELINE ###
#	Genome mapping with butter


############################    PARAMETERS    ##############################

# Input parameters
BASEDIR = "raw_compiled"
EXTENSION = ".txt.gz"

# Filtering parameters
FILTER_BASE = "A,T,C,G"
SIZE = "18,30"
TRIM = False
MIN_TRIM_LENGTH = 0

# Mapping parameters
TELO_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/telomere/bowtie/telomere"

# Other parameters
UTILS_DIR = "/nas/longleaf/home/sfrenk/scripts/util"

###############################################################################

# Get bowtie index for species
refs = {"elegans" : "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/bowtie/genome.fa", "remanei" : "/nas/longleaf/home/sfrenk/proj/seq/remanei/bowtie/genome.fa", "briggsae" : "/nas/longleaf/home/sfrenk/proj/seq/briggsae/WS263/bowtie/genome.fa"}

def get_ref(sample_name):

	species=re.search("elegans|briggsae|remanei", sample_name)
	if species is None:
		species="elegans"
	else:
		species = species.group(0)
	
	species_ref = refs[species]
	
	return(species_ref)

#if SPECIES not in refs:
#	print("ERROR: Unknown reference!")
#else:
#	REF = refs[SPECIES]

DATASET = re.search("([^/]*)/?$", BASEDIR).group(1)
SAMPLES = glob.glob(BASEDIR + "/*" + EXTENSION)
SAMPLES = [ re.search(BASEDIR + "/?([^/]+)" + EXTENSION, x).group(1) for x in SAMPLES ]


if len(SAMPLES) == 0:
	sys.exit("ERROR: no samples in base directory!")

rule all:
    input:
        expand("genome/{sample}_genome.bam", sample = SAMPLES),
        expand("results/{sample}.db", sample = SAMPLES),
        "results/" + DATASET + ".db"

rule filter_srna:
	input:
		BASEDIR + "/{sample}" + EXTENSION
	output:
		"filtered/{sample}.fa"
	params:
		filter_base = FILTER_BASE,
		size = SIZE,
		trim = TRIM,
		min_trim_length = MIN_TRIM_LENGTH,
		utils_dir = UTILS_DIR
	log:
		"logs/{sample}_filter.log"
	shell:
		"module add python; "
		"python3 {params.utils_dir}/small_rna_filter.py \
		-f {params.filter_base} \
		-s {params.size} \
		-t {params.trim} \
		-m {params.min_trim_length} \
		-o {output} \
		{input} > {log} 2>&1"

rule butter_mapping_genome:
	input:
		"filtered/{sample}.fa"
	output:
		bamfile = "genome/{sample}_genome.bam",
		bamidx = "genome/{sample}_genome.bam.bai"
	params:
		ref = lambda wildcards: get_ref('{sample}'.format(sample = wildcards.sample)),
		sample_name="{sample}"
	log:
		"logs/{sample}_map_genome.log"
	threads: 8
	shell:
		"module purge; " 
		"module add bowtie/1.1.2 perl samtools/0.1.19 && \
		butter --no_condense --aln_cores {threads} --max_rep 10000 --bam2wig none {input} {params.ref} > {log} 2>&1; "
		"mv filtered/{params.sample_name}.bam genome/{params.sample_name}_genome.bam; "
		"mv filtered/{params.sample_name}.bam.bai genome/{params.sample_name}_genome.bam.bai; "

rule map_telo_0_mismatch:
	input: "filtered/{sample}.fa"
	output:
		un = "fasta/{sample}_unmapped_0.fa",
		bam = "telo/{sample}_0.bam"
	params:
		telo_index = TELO_INDEX,
		al = "fasta/{sample}_0_al.fa",
		m = "fasta/{sample}_0_max.fa"
	threads: 8
	log:
		"logs/{sample}_map_0.log"
	shell:
		"module purge; " 
		"module add bowtie/1.1.2 samtools/1.8 && \
		bowtie -f --best --strata -M 1 -S -v 0 -p {threads} --al {params.al} --un {output.un} --max {params.m} {params.telo_index} {input} | samtools view -bh -F 4 - | samtools sort -o {output.bam} - > {log} 2>&1"

rule map_telo_1_mismatch:
	input: "fasta/{sample}_unmapped_0.fa"
	output:
		un = "fasta/{sample}_unmapped_1.fa",
		bam = "telo/{sample}_1.bam"
	params:
		telo_index = TELO_INDEX,
		al = "fasta/{sample}_1_al.fa",
		m = "fasta/{sample}_1_max.fa"
	threads: 8
	log:
		"logs/{sample}_map_1.log"
	shell:
		"module purge; " 
		"module add bowtie/1.1.2 samtools/1.8 && \
		bowtie -f --best --strata -M 1 -S -v 1 -p {threads} --al {params.al} --un {output.un} --max {params.m} {params.telo_index} {input} | samtools view -bh -F 4 - | samtools sort -o {output.bam} - > {log} 2>&1"

rule map_telo_2_mismatch:
	input: "fasta/{sample}_unmapped_1.fa"
	output:
		un = "fasta/{sample}_unmapped_2.fa",
		bam = "telo/{sample}_2.bam"
	params:
		telo_index = TELO_INDEX,
		al = "fasta/{sample}_2_al.fa",
		m = "fasta/{sample}_2_max.fa"
	threads: 8
	log:
		"logs/{sample}_map_2.log"
	shell:
		"module purge; " 
		"module add bowtie/1.1.2 samtools/1.8 && \
		bowtie -f --best --strata -M 1 -S -v 2 -p {threads} --al {params.al} --un {output.un} --max {params.m} {params.telo_index} {input} | samtools view -bh -F 4 - | samtools sort -o {output.bam} - > {log} 2>&1"

rule map_telo_3_mismatch:
	input: "fasta/{sample}_unmapped_2.fa"
	output:
		un = "fasta/{sample}_unmapped_3.fa",
		bam = "telo/{sample}_2.bam"
	params:
		telo_index = TELO_INDEX,
		al = "fasta/{sample}_3_al.fa",
		m = "fasta/{sample}_3_max.fa"
	threads: 8
	log:
		"logs/{sample}_map_3.log"
	shell:
		"module purge; " 
		"module add bowtie/1.1.2 samtools/1.8 && \
		bowtie -f --best --strata -M 1 -S -v 3 -p {threads} --al {params.al} --un {output.un} --max {params.m} {params.telo_index} {input} | samtools view -bh -F 4 - | samtools sort -o {output.bam} - > {log} 2>&1"

rule merge_telo_reads:
	input: "fasta/{sample}_unmapped_3.fa"
	output: "fasta/{sample}_all.fa"
	params:
		sample_name = "{sample}"
	run:
		shell('''if compgen -G "fasta/{params.sample_name}_[0-3]_*.fa" > /dev/null; then cat fasta/{params.sample_name}_[0-3]_*.fa > {output}; else touch {output}; fi''')

rule get_sample_info:
	input:
		telo_reads_file = "fasta/{sample}_all.fa",
		genome_bam = "genome/{sample}_genome.bam",
		bam_index = "genome/{sample}_genome.bam.bai"
	output:
		"sample_info/{sample}_info.txt"
	params:
		dataset = DATASET,
		sample_name = "{sample}"
	run:
		"module add samtools; "
		shell('''telo_reads=$(wc -l < {input.telo_reads_file}); module add samtools/1.8 && samtools view -c {input.genome_bam} | awk -v var="$telo_reads" '{{ print "{params.dataset}_{params.sample_name}\t{params.dataset}\t{params.sample_name}\t"$0"\t"var }}' > {output}''')

rule map_telo_reads_to_genome:
	input: 
		genome_bam = "genome/{sample}_genome.bam",
		fasta = "fasta/{sample}_all.fa"
	output: "telo/{sample}_genome.bam"
	params:
		idx = lambda wildcards: get_ref('{sample}'.format(sample = wildcards.sample))
	threads: 8
	shell:
		"module purge; " 
		"module add bowtie/1.1.2 samtools/1.8 && \
		if [[ $(wc -l < {input.fasta}) -eq 0 ]]; then samtools view -bH {input.genome_bam} > {output}; else bowtie -f --best --strata -a -S -v 0 -p {threads} {params.idx} {input.fasta} | samtools view -bh -F 4 - | samtools sort -o {output} -; fi"

rule convert_data:
	input: "telo/{sample}_genome.bam"
	output: 
		alignment_file = "results/{sample}_alignments.txt",
		reads_file = "results/{sample}_reads.txt"
	params:
		dataset = DATASET,
		sample_name = "{sample}",
		utils_dir = UTILS_DIR
	run:
		shell('''module add samtools/1.8 && if [[ $(samtools view {input} | wc -l) -eq 0 ]]; then touch {output.alignment_file}; else samtools view {input} | awk -v OFS="\t" -v d={params.dataset} -v s={params.sample_name} '{{print d"_"s"_"$1,$3,$4,1-($2/16),$6}}' > {output.alignment_file}; fi''')
		shell("module add python && python3 {params.utils_dir}/compile_data.py -d {params.dataset} -b fasta/{params.sample_name} -o {output.reads_file}")

rule store_results:
	input:
		alignment_file = "results/{sample}_alignments.txt",
		reads_file = "results/{sample}_reads.txt",
		sample_info = "sample_info/{sample}_info.txt"
	output:
		"results/{sample}.db"
	param:
		utils_dir = UTILS_DIR
	shell:
		"module add python; "
		"python3 {params.utils_dir}/ngs_store.py -d {output} -s {input.sample_info} -r {input.reads_file} -a {input.alignment_file}"

rule merge_db:
	input:
		expand("results/{sample}.db", sample = SAMPLES)
	output:
		"results/" + DATASET + ".db"
	params:
		utils_dir = UTILS_DIR
	shell:
		"module add python; "
		"python3 {parmas.utils_dir}/merge_db.py -o {output} {input}"
