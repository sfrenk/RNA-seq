#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("DESeq2"))

###############################################################################
# R script for simple ifferential gene expression analysis with DESeq2
###############################################################################

parser <- ArgumentParser(description = "basic differential expression with DESeq2. Input file should be a count table created using subread featureCounts. Note: the sample columns start at index 1")

parser$add_argument("file",
                    help = "input file")

parser$add_argument("-o", "--output",
                    help = "output file prefix (default: DESeq2_results)",
                    type = "character",
                    default = "DESeq2_results")

parser$add_argument("-c", "--control",
                    help = "column index (or indices) for control sample(s) as a comma-separated list",
                    type = "character")

parser$add_argument("-t", "--treatment",
                    help = "column index (or indices) for treatment sample(s) as a comma-separated list",
                    type = "character")

parser$add_argument("-f", "--total_counts",
                    help = "if normalizing to total library size (eg. for small RNA libraries) use this option followed by the total_mapped_reads.txt file",
                    default = "None")

parser$add_argument("-s", "--stringtie",
                    help = "count table was produced by stringtie",
                    action = "store_true")

args <- parser$parse_args()

count_file <- args$file

control_cols <- as.numeric(unlist(strsplit(args$control, ",")))
treatment_cols <- as.numeric(unlist(strsplit(args$treatment, ",")))

# experimental design
# 1. choose which columns to analyze. Note that the sample columns start at column number 6 from a featureCounts ouput file or column 2 from a merge_counts.py file
sample_cols <- c(control_cols, treatment_cols)

# 2. specify the condition for each sample (eg. "treatment" or "control")
conditions <- c(rep("control", length(control_cols)), rep("treatment", length(treatment_cols)))

################################################################################

if (args$stringtie) {
    
    # Stringtie output file
    rawcounts <- as.matrix(read.csv(count_file, row.names="transcript_id"))

} else {
    
    # Subread output file
    
    # Read in count data.
    rawcounts <- read.table(count_file, header = T, sep = "\t", row.names = 1)
    
    # Filter out rRNA reads
    rawcounts <- rawcounts[!grepl("rrn|rRNAinc", row.names(rawcounts)),]
    
    # Get rid of non-sample columns
    rawcounts <- rawcounts[,6:ncol(rawcounts)]
    
}


# Define experimental design:
rawcounts <- rawcounts[,sample_cols]

# Round counts (they may be decimal, eg if using --fraction option in subread)
rawcounts <- round(rawcounts, 0)

# Create sample_info object
names <- colnames(rawcounts)
sample_info <- data.frame(condition = conditions)
row.names(sample_info) <- names
print(sample_info)

# Quantifiy differential expression with DESeq
dds_count_table <- DESeqDataSetFromMatrix(countData = rawcounts,
                                          colData = sample_info,
                                          design = ~condition)
if (args$total_counts != "None"){
    
    # Create custom size factors
    totals <- read.table(args$total_counts)
    colnames(totals) <- c("sample", "count")
    mapped <- numeric()
    for (i in colnames(rawcounts)){
        
        # Get total number of mapped reads for each library
        mapped[i] <- totals[totals$sample == i, "count"]
    }
    
    # Get size factor for each library by dividing the total number of mapped reads by the average for the selected libraries.
    mapped_av <- mean(mapped)
    sizeFactors(dds_count_table) <- sapply(mapped, function(x) x/mapped_av)
}

dds <- DESeq(dds_count_table)
res <- results(dds, contrast=c("condition","treatment","control"))
res
res <- res[order(res$log2FoldChange),]
res <- na.omit(res)

# Output results
write.table(res, file = paste0(args$output, ".txt"), sep ="\t", col.names = NA, quote = F)

# MA plot
png(paste0(args$output, "_maPlot.png"))
plotMA(res)
dev.off()

# Dispersion plot
png(paste0(args$output, "_dispersionPlot.png"))
plotDispEsts(dds)
dev.off()

# PCA plot
rld <- rlog(dds, blind = FALSE)
png(paste0(args$output, "_PCAPlot.png"))
plotPCA(rld)
dev.off()

# Histogram of p values
png(paste0(args$output, "_pval_hist.png"))
hist(res$padj, breaks=100, col="blue", border="slateblue", main="p values")
dev.off()
