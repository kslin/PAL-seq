################################################################################
#DemultiplexAdjust.py
#Timothy J Eisen 2020 04 06
#This is a patch script for dealing with the fact that the intensity line isn't
#  in the same order as the fastq files. It iterates through each of the fastq 
#  files, asserting they are in the same order, and then iterates through the
#  intensity file, moving reads to the correct files. 
#  Code is now parallelized. 
################################################################################



import sys
from Bio import SeqIO
from collections import defaultdict as dd
import os
import pdb
import concurrent.futures
import gzip
sys.path.insert(1, '/lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqV3/tail-seq-scripts/')
import config

def genSet(FastqStrTuple):
	Fastq1Str, Fastq1Str = FastqStrTuple
	haFastq1 = open(Fastq1Str, 'r')
	haFastq2 = open(Fastq1Str, 'r')
	FullSet = set()
	genFastq1 = SeqIO.parse(haFastq1, "fastq") #generators
	genFastq2 = SeqIO.parse(haFastq2, "fastq")
	
	for record1 in genFastq1: #iterate thru fastq
		record2 = next(genFastq2) #same order
		FASTQ1_ID = config.fastq_header_to_ID(record1.id)
		FASTQ2_ID = config.fastq_header_to_ID(record2.id)
		assert FASTQ1_ID == FASTQ2_ID, 'FASTQ files out of order' #Are 
		FullSet.add(FASTQ1_ID)
	haFastq1.close()
	haFastq2.close()
	return(FullSet)


dictBarcodes = {}
IfileWrite = {}
dictSets = dd(set)
with open(sys.argv[1],'r') as haBarcodes:
	for line in haBarcodes: #creates a dict: (regex, filehandle)
		strBarcode    = line.split("\t")[0]
		strSampleName = line.split("\t")[1].strip()
		dictBarcodes[line.split("\t")[1].strip()] = (
			sys.argv[2] + '/' + strSampleName + '_1.txt', #just strings now.
			sys.argv[2] + '/' + strSampleName + '_2.txt')
		IfileWrite[line.split("\t")[1].strip()] = open(sys.argv[2] + '/' + strSampleName + '_HITS_MOD.txt', 'w+') #the only set of files that are written.


with concurrent.futures.ProcessPoolExecutor(max_workers = os.cpu_count()) as executor: #max workers is min(32, os.cpu_count())
 	for barcode, dataset in zip(dictBarcodes.keys(), executor.map(genSet, dictBarcodes.values())):
 		dictSets[barcode] = dataset

haIntensity = gzip.open(sys.argv[3], "rt") #open instensity, also gzipped. 
# haIntensity = open(sys.argv[3], 'r') #open instensity, not. 

for line in haIntensity:
	IlineID = config.intensity_line_to_ID(line.split())
	for identifier, dataset in dictSets.items():
		if IlineID in dataset:
			IfileWrite[identifier].write(line)
			break

#Close everything
for tups in IfileWrite.values():
	tups.close()
haIntensity.close()