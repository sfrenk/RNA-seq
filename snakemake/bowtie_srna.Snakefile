import glob
import re
import sys

############################    PARAMETERS    ##############################

# Input parameters
BASEDIR = ""
EXTENSION = ".fa.gz"

# Filtering parameters
FILTER_BASE = ["A", "T", "C", "G"]
SIZE = [22,22]
TRIM = False
MIN_TRIM_LENGTH = 0

# Mapping parameters
MULTI_FLAG = "-m 1"
MISMATCH = "-v 0"
REF = "genome"
REF_FASTA = "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/genome.fa"

# Counting parameters
ANTISENSE =  False
GTF = "/nas02/home/s/f/sfrenk/proj/seq/WS251/genes.gtf"

###############################################################################

SAMPLES = glob.glob(BASEDIR + "/*" + EXTENSION)
SAMPLES = [ re.search(BASEDIR + "/?([^/]+)" + EXTENSION, x).group(1) for x in SAMPLES ]

if len(SAMPLES) == 0:
	sys.exit("ERROR: no samples in base directory!")

rule all:
    input:
        "count/counts.txt"

rule filter_srna:
	input:
		BASEDIR + "/{sample}" + EXTENSION
	output:
		"filtered/{sample}.fa"
	params:
		filter_base = FILTER_BASE,
		size = SIZE,
		trim = TRIM,
		min_trim_length = MIN_TRIM_LENGTH
	log:
		"logs/{sample}_filter.log"
	shell:
		"small_rna_filter \
		-f {params.filter_base} \
		-s {params.size} \
		-t {params.trim} \
		-m {params.min_trim_length} \
		-o {output} \
		{input} > {log} 2>&1"

rule bowtie_index:
	input: REF_FASTA
	output: "index/genome.1.ebwt"
	params:
		output_name = "index/genome"
	shell:
  		"bowtie-build {input} {params.output_name}"

rule bowtie_mapping:
	input:
		reads_file = "filtered/{sample}.fa",
		idx_file = "index/genome.1.ebwt"
	output: "bowtie_out/{sample}.sam"
	params:
		idx = "index/genome",
		multi_flag = MULTI_FLAG,
		mismatch = MISMATCH
	log:
		"logs/{sample}_map.log"
	threads: 4
	shell: "bowtie -f --best --strata -S \
		-p {threads} \
		{params.multi_flag} \
		-v {params.mismatch} \
		{params.idx} \
		{input.reads_file} {output} > {log} 2>&1" 

rule convert_to_bam:
	input: "bowtie_out/{sample}.sam"
	output: "bam/{sample}.bam"
	shell:
		"samtools view -bh -F 4 {input} | samtools sort -o {output} -"

rule index_bam:
	input: "bam/{sample}.bam"
	shell:
		"samtools index {input}"

if REF != "genome":
	rule count_reads:
		input:
			bam_file = "bam/{sample}.bam"
		output:
			"count/{sample}_counts.txt"
		threads: 4
		run:
			if ANTISENSE:
				shell('''samtools view -f 16 {input.bam_file} | awk -F '\t' '$2 == "16" ' | grep -v @ | cut -f 3 | sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' > {output}''')
			else:
				shell('''samtools view {input.bam_file} | awk -F '\t' '$2 == "16" ' | grep -v @ | cut -f 3 | sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' > {output}''')

	rule merge_counts:
		input:
			expand("count/{sample}.txt", sample = SAMPLES)
		output:
			"count/counts.txt"
		threads: 4
		run:
			"python /proj/ahmedlab/steve/seq/util/merge_counts.py -o {output} ./count"

else:
	# Count all samples at the same time with subread
	rule count_reads_genome:
		input:
			expand("bam/{sample}.bam", sample = SAMPLES)
		output:
			"count/counts.txt"
		params:
			gtf = GTF
		log:
			"logs/count.log"
		threads: 4
		run:
			if ANTISENSE:
				shell("featureCounts -M -s 2 -a {params.gtf} -o {output} -T {threads} -t exon -g gene_name {input} > {log} 2>&1")
			else:
				shell("featureCounts -M -a {params.gtf} -o {output} -T {threads} -t exon -g gene_name {input} > {log} 2>&1")
