import glob
import re
import sys

############################    PARAMETERS    ##############################

# Input parameters
BASEDIR = ""
EXTENSION = ""
PAIRED = False

# Trimming parameters
ADAPTERS="~/proj/seq/bbmap/adapters.fa"

# Mapping parameters
INDEX = "/nas02/home/s/f/sfrenk/proj/seq/WS251/genome/hisat2/genome"
REMOVE_RDNA = False # Set to True if using ribosome-depleted samples
RDNA_BED = "/nas02/home/s/f/sfrenk/proj/seq/rdna/rdna_loci.bed"


# Counting parameters
COUNT_METHOD="subread"
GTF = "/nas02/home/s/f/sfrenk/proj/seq/WS251/genes.gtf"

###############################################################################

SAMPLE_FILES = glob.glob(BASEDIR + "/*" + EXTENSION)
SAMPLES = [ re.search(BASEDIR + "/?([^/]+)" + EXTENSION, x).group(1) for x in SAMPLE_FILES ]

# Paired end files must end in _1/_2 where 1 and 2 denote forward and reverse reads respectively. 
if PAIRED:
	SAMPLES = list(set([re.search("(^.+)_[12]$", x).group(1) for x in SAMPLES]))
	#SAMPLES = list(set([re.search("(^.+)_[rR]?[12]$", x).group(1) for x in SAMPLES]))

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


#### Mapping ####

if PAIRED:

	rule fastqc_pre_trim:
		input:
			read1 = BASEDIR + "/{sample}_1" + EXTENSION,
			read2 = BASEDIR + "/{sample}_2" + EXTENSION
		output:
			html1 = "metrics/fastq/{sample}_1_fastqc.html",
			html2 = "metrics/fastq/{sample}_2_fastqc.html",
			zip1 = "metrics/fastq/{sample}_1_fastqc.zip",
			zip2 = "metrics/fastq/{sample}_2_fastqc.zip"
		threads: 4
		log:
			"logs/{sample}_fastqc_pre_trim.log"
		shell:
			"module add fastqc; "
			"fastqc -o metrics/fastq {input.read1} {input.read2} &> {log}"

	rule trim:
		input:
			read1 = BASEDIR + "/{sample}_1" + EXTENSION,
			read2 = BASEDIR + "/{sample}_2" + EXTENSION
		output:
			out1 = "trimmed/{sample}_1.fastq.gz",
			out2 = "trimmed/{sample}_2.fastq.gz"
		params:
			adapter_file = ADAPTERS
		threads: 1
		log:
			"logs/{sample}_trim.log"
		shell:
			"module add bbmap; "
			"bbduk.sh -Xmx4g in1={input.read1} in2={input.read2} out1={output.out1} out2={output.out2} ref={params.adapter_file} ktrim=r overwrite=true k=23 maq=20 mink=11 hdist=1 &> {log}"

	rule fastqc_post_trim:
		input:
			read1 = "trimmed/{sample}_1_trimmed.fastq",
			read2 = "trimmed/{sample}_2_trimmed.fastq"
		output:
			html1 = "metrics/fastq/{sample}_1_trimmed_fastqc.html",
			html2 = "metrics/fastq/{sample}_2_trimmed_fastqc.html",
			zip1 = "metrics/fastq/{sample}_1_trimmed_fastqc.zip",
			zip2 = "metrics/fastq/{sample}_2_trimmed_fastqc.zip"
		threads: 4
		log:
			"logs/{sample}_fastqc_post_trim.log"
		shell:
			"module add fastqc; "
			"fastqc -o metrics/fastq {input.read1} {input.read2} &> {log}"

	rule hisat2_mapping:
		input:
			trimmed1 = "trimmed/{sample}_1.fastq.gz",
			trimmed2 = "trimmed/{sample}_2.fastq.gz"
		output:
			"hisat2_out/{sample}.sam"
		params:
			idx_base = INDEX,
		log:
			"logs/{sample}_map.log"
		threads: 8
		shell: 
			"module add hisat2; "
			"hisat2 --max-intronlen 12000 --dta --no-mixed --no-discordant -p {threads} -x {params.idx_base} -1 {input.trimmed1} -2 {input.trimmed2} -S {output} &> {log}"


else:

	rule fastqc_pre_trim:
		input:
			BASEDIR + "/{sample}" + EXTENSION
		output:
			html = "metrics/fastq/{sample}_fastqc.html",
			zipfile = "metrics/fastq/{sample}_fastqc.zip"
		threads: 4
		log:
			"logs/{sample}_fastqc_pre_trim.log"
		shell:
			"module add fastqc; "
			"fastqc -o metrics/fastq {input} &> {log}"

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
			"module add bbmap; "
			"bbduk.sh -Xmx4g in={input} out={output} ref={params.adapter_file} ktrim=r overwrite=true k=23 maq=20 mink=11 hdist=1 &> {log}"
	
	rule fastqc_post_trim:
		input:
			"trimmed/{sample}_trimmed.fastq"
		output:
			html = "metrics/fastq/{sample}_trimmed_fastqc.html",
			zipfile = "metrics/fastq/{sample}_trimmed_fastqc.zip"
		threads: 4
		log:
			"logs/{sample}_fastqc_post_trim.log"
		shell:
			"module add fastqc; "
			"fastqc -o metrics/fastq {input} &> {log}"

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
		shell:
			"module add hisat2; " 
			"hisat2 --max-intronlen 12000 --dta --no-mixed --no-discordant -p {threads} -x {params.idx_base} -U {input} -S {output} &> {log}"


#### rDNA removal ####

if REMOVE_RDNA == True:

	rule convert_to_bam:
		input:
			"hisat2_out/{sample}.sam"
		output:
			"temp/{sample}.temp.bam"
		params:
			rdna_bed = RDNA_BED
		shell:
			"module add samtools; "
			"samtools view -bh -L {params.rdna_bed} {input} -U {output}"


	rule sort_bam:
		input:
			"temp/{sample}.temp.bam"
		output:
			"bam/{sample}.bam"
		params:
			rdna_bed = RDNA_BED
		shell:
			"module add samtools; "
			"samtools view -bh -F 4 {input} | samtools sort -o {output} -"

else:

	rule convert_to_bam:
		input:
			"hisat2_out/{sample}.sam"
		output:
			"bam/{sample}.bam"
		shell:
			"module add samtools; "
			"samtools view -bh -F 4 {input} | samtools sort -o {output} -"

rule index_bam:
	input:
		"bam/{sample}.bam"
	output:
		"bam/{sample}.bam.bai" 
	shell:
		"module add samtools; "
		"samtools index {input}"

rule bam_metrics:
	input:
		bamfile = "bam/{sample}.bam",
		bamidx = "bam/{sample}.bam.bai"
	output:
		"metrics/bam/{sample}/qualimapReport.html"
	params:
		output_dir="metrics/bam/{sample}"
	log:
		"logs/{sample}_bam_metrics.log"
	shell:
		"module add qualimap; "
		"qualimap bamqc -bam {input.bamfile} -outdir {params.output_dir} -gd HUMAN &> {log}"


#### QC compilation ####

if PAIRED:
	rule multiqc:
		input:
			fastq_pre_trim_1 = expand("metrics/fastq/{sample}_1_fastqc.html", sample = SAMPLES),
			fastq_pre_trim_2 = expand("metrics/fastq/{sample}_2_fastqc.html", sample = SAMPLES),
			bam = expand("metrics/bam/{sample}/qualimapReport.html", sample = SAMPLES),
			fastq_post_trim_1 = expand("metrics/fastq/{sample}_1_trimmed_fastqc.html", sample = SAMPLES),
			fastq_post_trim_2 = expand("metrics/fastq/{sample}_2_trimmed_fastqc.html", sample = SAMPLES),
			picard = expand("metrics/picard/{sample}_metrics.txt", sample = SAMPLES)
		output:
			"multiqc_report.html"
		log:
			"logs/multiqc.log"
		shell:
			"module purge; "
			"module add multiqc; "
			"multiqc -f metrics &> {log}"
else:
	rule multiqc:
		input:
			fastq_pre_trim = expand("metrics/fastq/{sample}_fastqc.html", sample = SAMPLES),
			fastq_post_trim = expand("metrics/fastq/{sample}_trimmed_fastqc.html", sample = SAMPLES),
			bam = expand("metrics/bam/{sample}/qualimapReport.html", sample = SAMPLES),
			picard = expand("metrics/picard/{sample}_metrics.txt", sample = SAMPLES)
		output:
			"multiqc_report.html"
		log:
			"logs/multiqc.log"
		shell:
			"module purge; "
			"module add multiqc; "
			"multiqc -f metrics &> {log}"


#### Assigning reads to features ####

if COUNT_METHOD == "subread":
	rule count:
		input:
			bamfiles = expand("bam/{sample}.bam", sample = SAMPLES),
			bamfile_idx = expand("bam/{sample}.bam.bai", sample = SAMPLES),
			gtf = GTF #gtf = "stringtie/stringtie_merged.gtf"
		output:
			"count/counts.txt"
		log:
			"logs/count.log"
		threads: 8
		shell:
			"module add subread; "
			"featureCounts -a {input.gtf} -M --fraction -o {output} -T {threads} -t exon -g gene_name {input.bamfiles} &> {log}"

else:
	rule stringtie_make:
		input:
			bamfile = "bam/{sample}.bam",
			bamidx = "bam/{sample}.bam.bai"
		output:
			"stringtie/{sample}.gtf"
		params:
			gtf = GTF
		log:
			"logs/{sample}_stringtie_make.log"
		threads: 8
		shell:
			"module add stringtie; "
			"stringtie {input.bamfile} -p {threads} -o {output} -G {params.gtf} &> {log}"

	rule stringtie_merge:
		input:
			expand("stringtie/{sample}.gtf", sample = SAMPLES)
		output:
			"stringtie/stringtie_merged.gtf"
		params:
			gtf = GTF
		log:
			"logs/stringtie_merge.log"
		threads: 8
		run:
			"module add stringtie; "
			"ls stringtie/*.gtf > gtf_files.txt; "
			"stringtie --merge -p {threads} -G {params.gtf} -o {output} gtf_files.txt &> {log}"

	rule count:
		input:
			bamfile = "bam/{sample}.bam",
			gtf = "stringtie/stringtie_merged.gtf"
		output:
			"stringtie_count/{sample}/e2t.ctab"
		log:
			"logs/{sample}_stringtie_count.log"
		threads: 8
		shell:
			"module add stringtie; "
			"stringtie {input.bamfile} -p {threads} -B -e -o {output} -G {input.gtf} &> {log}"
