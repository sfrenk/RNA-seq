import glob
import re
import sys

############################    PARAMETERS    ##############################

# Input parameters
BASEDIR = ""
EXTENSION = ".fa.gz"

# Filtering parameters
FILTER_BASE = "A,T,C,G"
SIZE = "18,30"
TRIM = False
MIN_TRIM_LENGTH = 0

# Mapping parameters
MULTI_FLAG = "-M 1"
MISMATCH_FLAG = "-v 0"
REF = "genome"

# Counting parameters
ANTISENSE =  False
GTF = "/nas02/home/s/f/sfrenk/proj/seq/WS251/genes.gtf"

###############################################################################

# Get bowtie index
if REF == "genome":
	BOWTIE_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/bowtie/genome"
elif REF == "transposons":
	BOWTIE_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/transposons/bowtie/transposon"
elif REF == "mrna":
	BOWTIE_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/WS251/mrna/bowtie/mrna"
elif REF == "mirna":
	BOWTIE_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/mirna/bowtie/mirna"
elif REF == "repeats":
	BOWTIE_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/repeats/repbase/simple_repeats/bowtie/repeats"
else:
	print("ERROR: Unknown reference!")

# Get samples
SAMPLES = glob.glob(BASEDIR + "/*" + EXTENSION)
SAMPLES = [ re.search(BASEDIR + "/?([^/]+)" + EXTENSION, x).group(1) for x in SAMPLES ]

# Pipeline
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

rule bowtie_mapping:
	input:
		"filtered/{sample}.fa"
	output: "bowtie_out/{sample}.sam"
	params:
		idx = BOWTIE_INDEX,
		multi_flag = MULTI_FLAG,
		mismatch_flag = MISMATCH_FLAG
	log:
		"logs/{sample}_map.log"
	threads: 8
	shell:
		"module purge; " 
		"module add bowtie/1.1.2; "
		"bowtie -f --best --strata -S \
		-p {threads} \
		{params.multi_flag} \
		{params.mismatch_flag} \
		{params.idx} \
		{input} {output} > {log} 2>&1" 

rule convert_to_bam:
	input: "bowtie_out/{sample}.sam"
	output: "bam/{sample}.bam"
	shell:
		"samtools view -bh -F 4 {input} | samtools sort -o {output} -"

rule index_bam:
	input:
		"bam/{sample}.bam"
	output:
		"bam/{sample}.bam.bai"
	shell:
		"samtools index {input}"

if REF != "genome":
	rule count_reads:
		input:
			bamfile = "bam/{sample}.bam",
			bamidx = "bam/{sample}.bam.bai"
		output:
			"count/{sample}.txt"
		threads: 1
		run:
			if ANTISENSE:
				shell("samtools view -f 16 {input.bamfile} | cut -f 3 | sort | uniq -c |  sed 's/^[ \t]*//g' | awk -v OFS='\t' '{{print $1,$2}}' > {output}")
			else:
				shell("samtools view {input.bamfile} | cut -f 3 | sort | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS='\t' '{{print $1,$2}}' > {output}")

	rule merge_counts:
		input:
			expand("count/{sample}.txt", sample = SAMPLES)
		output:
			"count/counts.txt"
		threads: 1
		shell:
			"python3 /proj/ahmedlab/steve/seq/util/merge_counts.py -o {output} {input}"

else:
	# Count all samples at the same time with subread
	rule count_reads_genome:
		input:
			bamfiles = expand("bam/{sample}.bam", sample = SAMPLES),
			bamidx = expand("bam/{sample}.bam.bai", sample = SAMPLES)
		output:
			"count/counts.txt"
		params:
			gtf = GTF
		log:
			"logs/count.log"
		threads: 4
		run:
			if ANTISENSE:
				shell("featureCounts -M -s 2 -a {params.gtf} -o {output} -T {threads} -t exon -g gene_name {input.bamfiles} > {log} 2>&1")
			else:
				shell("featureCounts -M -a {params.gtf} -o {output} -T {threads} -t exon -g gene_name {input.bamfiles} > {log} 2>&1")
