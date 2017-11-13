import glob
import re
import sys

COUNT_METHOD="subread"
BASEDIR = "fastq"

GTF = "/nas02/home/s/f/sfrenk/proj/seq/WS251/genes.gtf"
INDEX = "/nas02/home/s/f/sfrenk/proj/seq/WS251/genome/hisat2/genome"
ADAPTERS="/nas/longleaf/apps/bbmap/37.62/bbmap/resources/adapters.fa"
EXTENSION = ""

SAMPLES = glob.glob(BASEDIR + "/*" + EXTENSION)
SAMPLES = [ re.search(BASEDIR + "/?([^/]+)" + EXTENSION, x).group(1) for x in SAMPLES ]

if len(SAMPLES) == 0:
	sys.exit("ERROR: no samples in base directory!")

if COUNT_METHOD == "subread":
	rule all:
		input:
			"count/counts.txt"
elif COUNT_METHOD == "stringtie":
	rule all:
		input:
			expand("stringtie_count/{sample}/e2t.ctab", sample = SAMPLES)
else:
	sys.exit("ERROR: Invalid COUNT_METHOD option. Choose subread or stringtie.")

rule trim:
	input:
		BASEDIR + "/{sample}" + EXTENSION
	output:
		"trimmed/{sample}.fastq.gz"
	params:
		adapter_file = ADAPTERS
	threads: 1
	log:
		"logs/{sample}_trim.log"
	shell:
		"bbduk.sh -Xmx4g in={input} out={output} ref={params.adapter_file} ktrim=r overwrite=true k=23 maq=20 mink=11 hdist=1 > {log} 2>&1"

rule hisat2_mapping:
	input:
		"trimmed/{sample}.fastq.gz"
	output:
		"hisat2_out/{sample}.sam"
	params:
		idx_base = INDEX,
	log:
		"logs/{sample}_map.log"
	threads: 8
	shell: "hisat2 --max-intronlen 12000 --dta --no-mixed --no-discordant -p {threads} -x {params.idx_base} -U {input} -S {output} > {log} 2>&1"

rule convert_to_bam:
	input:
		"hisat2_out/{sample}.sam"
	output:
		"bam/{sample}.bam"
	shell:
		"samtools view -bh -F 4 {input} | samtools sort -o {output} -"

rule index_bam:
	input: "bam/{sample}.bam"
	shell:
		"samtools index {input}"

if COUNT_METHOD == "subread":
	rule count:
		input:
			bamfiles = expand("bam/{sample}.bam", sample = SAMPLES),
			gtf = GTF #gtf = "stringtie/stringtie_merged.gtf"
		output:
			"count/counts.txt"
		log:
			"logs/count.log"
		threads: 8
		shell:
			"featureCounts -a {input.gtf} -M --fraction -o {output} -T {threads} -t exon -g gene_name {input.bamfiles} > {log} 2>&1"

else:
	rule stringtie_make:
		input:
			"bam/{sample}.bam"
		output:
			"stringtie/{sample}.gtf"
		params:
			gtf = GTF
		threads: 8
		shell:
			"stringtie {input} -p {threads} -o {output} -G {params.gtf}"

	rule stringtie_merge:
		input:
			expand("stringtie/{sample}.gtf", sample = SAMPLES)
		output:
			"stringtie/stringtie_merged.gtf"
		params:
			gtf = GTF
		threads: 8
		run:
			shell("ls stringtie/*.gtf > gtf_files.txt")
			shell("stringtie --merge -p {threads} -G {params.gtf} -o {output} gtf_files.txt")

	rule count:
		input:
			bamfile = "bam/{sample}.bam",
			gtf = "stringtie/stringtie_merged.gtf"
		output:
			"stringtie_count/{sample}/e2t.ctab"
		log:
			"logs/{sample}stringtie_count.log"
		threads: 8
		shell:
			"stringtie {input.bamfile} -p {threads} -B -e -o {output} -G {input.gtf}> {log} 2>&1"
