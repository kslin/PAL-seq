
#####################################################
#  Demultiplexing script for PALseq V3 pipeline
#  Expects a six nucleotide barcode as the first 6 nucleotides
#  of FASTQ1. Does not check chastity filter. Allows 1 error.
#  Second input is FASTQ2.
#  Third input is a tab-delimited list of barcodes and sample names
#  Uses Biopython
#  Timothy J Eisen
#  2019 09 19
#####################################################

import gzip
from Bio import SeqIO
import sys
import regex

#Generate the barcode dictionary
dictBarcodes = {}
with open(sys.argv[1],'r') as haBarcodes:
	for line in haBarcodes: #creates a dict: (regex, filehandle)
		strBarcode    = line.split("\t")[0]
		strSampleName = line.split("\t")[1].strip()
		dictBarcodes[line.split("\t")[1].strip()] = \
			(regex.compile('(%s){e<=1}' % (strBarcode)), #gen regex
			open(strSampleName + '_1.txt', 'w+'),
			open(strSampleName + '_2.txt', 'w+')) #open the file handles.

#iterator dictionary
genBarcodes = {}
haFastq1 = gzip.open(sys.argv[2], "rt") #open fastq1
haFastq2 = gzip.open(sys.argv[3], "rt") #open fastq2
genFastq1 = SeqIO.parse(haFastq1, "fastq") #generators
genFastq2 = SeqIO.parse(haFastq2, "fastq")

for record1 in genFastq1: #iterate thru fastq
	record2 = next(genFastq2) #same order
	for strSampleName, tupVal in dictBarcodes.items(): 
		if tupVal[0].match(str(record1.seq[:6])) is not None:
			tupVal[1].write(record1.format("fastq"))
			tupVal[2].write(record2.format("fastq"))
			break

#close all files.
haFastq1.close()
haFastq2.close()
[tups[1].close() for tups in dictBarcodes.values()]
[tups[2].close() for tups in dictBarcodes.values()]
