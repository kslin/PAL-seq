
#####################################################
#  Demultiplexing script for PALseq V3 pipeline
#  Expects a six nucleotide barcode as the first 6 nucleotides
#  of FASTQ1. Does not check chastity filter. Allows 1 error.
#  First input is barcode file, barcode tab sample name
#  Second input is FASTQ1.
#  Third input is FASTQ2.
#  Fourth input is the intesity file. 
#  Fifth output directory
#  Uses Biopython
#  Runtime of ~2.2 hrs on a lane of a rapid run
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
			open(sys.argv[5] + '/' + strSampleName + '_1.txt', 'w+'),
			open(sys.argv[5] + '/' + strSampleName + '_2.txt', 'w+'),
			open(sys.argv[5] + '/' + strSampleName + '_hits.txt', 'w+')) #open the file handles.

#iterator dictionary
genBarcodes = {}
haFastq1 = gzip.open(sys.argv[2], "rt") #open fastq1
haFastq2 = gzip.open(sys.argv[3], "rt") #open fastq2
haIntensity = open(sys.argv[4], "rt") #open fastq2
genFastq1 = SeqIO.parse(haFastq1, "fastq") #generators
genFastq2 = SeqIO.parse(haFastq2, "fastq")
genIntensity = (i for i in haIntensity)
for record1 in genFastq1: #iterate thru fastq
	record2 = next(genFastq2) #same order
	recordIntensity = next(genIntensity)
	for strSampleName, tupVal in dictBarcodes.items(): 
		if tupVal[0].match(str(record1.seq[:6])) is not None:
			tupVal[1].write(record1.format("fastq"))
			tupVal[2].write(record2.format("fastq"))
			tupVal[3].write(recordIntensity)
			break

#close all files.
haFastq1.close()
haFastq2.close()
haIntensity.close()
for tups in dictBarcodes.values():
	tups[1].close()
	tups[2].close()
	tups[3].close()
