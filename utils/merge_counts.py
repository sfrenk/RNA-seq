#!/usr/bin/env python

###############################################################################
# Use this script to merge all count.txt files in a count directory into a count table. The resulting table can be loaded straight into DESeq2 for analysis
###############################################################################

import pandas as pd
import os
import argparse

# Get target directory

parser = argparse.ArgumentParser(description = "Merge multiple count files into a single count table")

parser.add_argument("count_files", help = "Count files to merge", nargs = "+")
parser.add_argument("-o", "--output", help = "output file name (default: count_table.txt", default = "count_table.txt")

args = parser.parse_args()

# Read in first count file. All further counts will be appended to this data frame

merged = pd.read_csv(args.count_files[0], sep="\t", header=None, names=[args.count_files[0].replace("_counts.txt", ""), 'gene'])
merged = merged[['gene', args.count_files[0].replace("_counts.txt", "")]]

# Remove the first count file from the list so that it doesn't get counted twice

count_files = args.count_files[1:]

# Merge the other count files to the first one

for i in count_files:
	count_file = pd.read_csv(i, sep="\t", header=None, names=[i.replace("_counts.txt", ""), 'gene'])
	merged = merged.merge(count_file, how='left', on='gene')

# Get rid of the row containing unmapped reads

merged = merged[merged.gene != '*']

# If a particular transcript was detected in some samples but not others, the samples in which the transcript was absent will report an "NA" value after merging of the count files. This needs to be changed to "0" for downstream analysis

merged = merged.fillna(0)

# Save the count table

merged.to_csv(args.output, sep='\t', index=False)
