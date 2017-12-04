import glob
import re

REF = "genome"
ANTISENSE = True
BAMDIR = "bam"
GTF = "/nas02/home/s/f/sfrenk/proj/seq/WS251/genes.gtf"

SAMPLES = glob.glob(BAMDIR + "/*.bam")
SAMPLES = [ re.search(BAMDIR + "/?([^/]+).bam", x).group(1) for x in SAMPLES ]

if REF != "genome":
	rule count_reads:
		input:
			"bam/{sample}.bam"
		output:
			"count/{sample}.txt"
		threads: 1
		run:
			if ANTISENSE:
				shell("samtools view -f 16 {input} | cut -f 3 | sort | uniq -c |  sed 's/^[ \t]*//g' | awk -v OFS='\t' '{print $1,$2}' > {output}")
			else:
				shell("samtools view {input} | cut -f 3 | sort | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS='\t' '{print $1,$2}' > {output}")

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