#!/usr/bin/env python3

###############################################################################
# small_rna_filter.py version 3.0
#
#
# Filter reads from a fastq, fasta or raw file and output in any format.
# Optionally filter based on sequence composition (eg. to extract 22Gs)
###############################################################################

import os
import argparse
import gzip
import re
import sys

def determine_file_type(file):

	'''
	Determine file type based on extension
	'''

	mode = ""

	if re.search("\.fastq\.gz$|\.fastq$|\.fq$|\.fq\.gz$", file):
		mode = "fastq"
	elif re.search("\.fasta\.gz$|\.fasta$|\.fa$|\.fa\.gz$", file):
		mode = "fasta"
	elif re.search("\.txt\.gz$|\.txt$", file):
		mode = "txt"
	else:
		sys.exit("ERROR: Unrecognized file type: " + str(file))

	return(mode)


def trim_3prime(seq, base, min_trim_length = 0):

	'''
	Trim any 3' <base> nucleotides

	Arguments:

	seq -- Sequence string

	Returns:

	seq -- trimmed sequence
	count -- number of nucleotides removed
	'''
	count = 0
	while seq.upper().endswith(base.upper()) and len(seq) > min_trim_length:
		# Sequentially remove all 3' A nucleotides
		seq = seq[:-1]
		count += 1

	return(seq, count)

def trim_3prime_single(seq, base):
	'''
	If read ends with base, remove the last nucleotide and return the read. Otherwise, return nothing
	'''
	if seq.upper().endswith(base.upper()):
		seq = seq[:-1]
		return(seq)
	else:
		return(None)

def write_output(output_file, output_mode, header = None, seq = None, line3 = "+", qual = None, counter = ""):

	'''
	Write read to file in specified format
	'''

	if output_mode == "fastq":
		if header == None:
			header = "@Seq_" + str(counter)
		elif header[0] != "@":
			header = "@" + header

		if qual == None:
			qual = "I" * len(seq)
		else:
			# Trim the quality to match any trimming that happened to the sequence
			qual = qual[:len(seq)]

		output_file.write(header + "\n" + seq + "\n" + line3 + "\n" + qual + "\n")

	elif output_mode == "fasta":
		if header == None:
			header = ">Seq_" + str(counter)
		elif header[0] != ">":
			header = ">" + header

		output_file.write(header + "\n" + seq + "\n")

	elif output_mode == "txt":

		output_file.write(seq + "\n")


def filter_test(seq, filter_params):

	'''
	Return True if sequence fulfills filtering requirements, otherwise False
	'''

	if seq and len(seq) <= filter_params["max_size"] and len(seq) >= filter_params["min_size"]:
		if seq[0].upper() in filter_params["first_base"]:
			return True

	return False

def parse_input(input_file, output_file, input_mode, output_mode, filter_params, trim = None, trim_single = None, min_trim_length = 0):
	input_counter = 0
	output_counter = 0
	trim_counter = 0

	if input_mode == "fastq":

		while True:
			
			try:
				# Note: have to remove "@" from the beginning of fastq header 
				header = next(input_file).strip()[1:]
				seq = next(input_file).strip()

				if trim is not None:
					seq, count = trim_3prime(seq, base = trim, min_trim_length = min_trim_length)
					trim_counter += count

				if trim_single is not None:
					seq = trim_3prime_single(seq, trim_single)

				line3 = next(input_file).strip()
				qual = next(input_file).strip()

				if filter_test(seq, filter_params):

					output_counter += 1
					write_output(output_file, output_mode, header, seq, line3, qual, counter = output_counter)
					
				input_counter += 1
			
			except StopIteration:
				break


	elif input_mode == "fasta":
		
		while True:
			
			try:
				header = next(input_file).strip()
				seq = next(input_file).strip()
				
				if trim is not None:
					seq, count = trim_3prime(seq, base = trim, min_trim_length = min_trim_length)
					trim_counter += count

				if trim_single is not None:
					seq = trim_3prime_single(seq, trim_single)

				if filter_test(seq, filter_params):
					output_counter += 1
					write_output(output_file, output_mode, header, seq, counter = output_counter)
					
				input_counter += 1

			except StopIteration:
				break

	elif input_mode == "txt":

		while True:

			try:
				seq = next(input_file).strip()

				if trim is not None:
					seq, count = trim_3prime(seq, base = trim, min_trim_length = min_trim_length)
					trim_counter += count

				if trim_single is not None:
					seq = trim_3prime_single(seq, trim_single)

				if filter_test(seq, filter_params):
					output_counter += 1
					write_output(output_file, output_mode, seq = seq, counter = output_counter)

				input_counter += 1

			except StopIteration:
				break

	return(input_counter, output_counter, trim_counter)


def __main__():

	parser = argparse.ArgumentParser(description="Extract 22g and/or 21u rnas from fastq, fasta or raw text file")

	parser.add_argument("infile", help = "Input file. File type will be detected based on the first line of the file. File may be in gzipped format.")
	parser.add_argument("-o", "--output", help = "Output file.", required = True)
	parser.add_argument("-f", "--filter", help = "Filter reads based on the first nucleotide (eg. to select 22g RNAs, use options -f g -s 22,22", default = "A,T,G,C")
	parser.add_argument("-s", "--size", help = "Specify size range of reads to keep by providing min and max length seperated by a comma (eg. to keep reads between 19 and 24 nucleotides (inclusive), use '-s 19,24' (Default size range: 18 to 30 nucleotides)", default = "18,30")
	parser.add_argument("-t", "--trim", help = "trim 3' nucleotides of this base (eg. -t A to remove any 3' A nucleotides)", default = None)
	parser.add_argument("-m", "--min_trim_length", help = "Stop trimming 3' nucleotides when read has this many nucleotides", default = 0, type = int)
	parser.add_argument("-u", "--trim_single", help = "Only return reads ending in this base, with the final nucleotide removed (eg. remove the final 'T' to get potential CSR-1 RNAs)", default = None)
	

	args = parser.parse_args()

	# Define starting nucleotide

	args.filter = args.filter.upper()
	first_base = args.filter.split(",")

	# Define minimum and maximum read size

	min_size = int(args.size.split(",")[0])
	max_size = int(args.size.split(",")[1])

	# Compile filter parameters

	filter_params = { "first_base" : first_base, "min_size" : min_size, "max_size" : max_size }

	# Determine input/output file types

	input_mode = determine_file_type(args.infile)
	output_mode = determine_file_type(args.output)

	print("##### small_rna_filter.py #####")

	if args.trim_single is not None:
		print("Single trim mode: Returning reads ending in " + args.trim_single + ", with the final nucleotide removed.")

	# Input file may or may not be gzipped

	if args.infile.endswith(".gz"):
		input_file = gzip.open(args.infile, mode = "rt")
	else:
		input_file = open(args.infile, "r")

	# Perform filtering

	output_file = open(args.output, "w")
	
	input_count, output_count, trim_count = parse_input(input_file, output_file, input_mode, output_mode, filter_params, trim = args.trim, trim_single = args.trim_single, min_trim_length = args.min_trim_length)

	input_file.close()
	output_file.close()

	print("Input reads: " + str(input_count))
	if args.trim is not None:
		print("Removed " + str(trim_count) + " 3' " + args.trim + " nucleotides")
	print("Output reads: " + str(output_count))

	print("#####	Finished	#####")

if __name__ == "__main__":
	__main__()
