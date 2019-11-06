#Opens the fastq files from their original positions on solexa public as tarfiles, then uses Biopython to parse
#Must be run with python 2.7
#python2.7 tail-seq-scripts/dephase.py /lab/solexa_public/Bartel/160311_WIGTC-HISEQA_C73V2ACXX/QualityScore/CTAGTC-s_4_1_sequence.txt.tar.gz
#Include extra packages in requirements file??

from Bio import SeqIO, bgzf
import sys
import tarfile
import re
from gzip import open as gzopen


fastq1TAR=tarfile.open(name=sys.argv[1], mode='r:gz')
fastq1open=fastq1TAR.extractfile(fastq1TAR.next())

fastq2TAR=tarfile.open(name=sys.argv[2], mode='r:gz')
fastq2open=fastq2TAR.extractfile(fastq2TAR.next())

counter = 0

dictPhase = {'objPhase1':re.compile("A..."),'objPhase2':re.compile("GA.."),'objPhase3':re.compile("CGA."),'objPhase4':re.compile("TCGA")}
dictTrimPhase = {'objPhase1':9,'objPhase2':10,'objPhase3':11,'objPhase4':12}

def trim_adaptors(records, dictPhase, dictTrimPhase):
	"""Trims adapter sequencing depending on phasing.

	This is a generator function, the records argument should
	be a list or iterator returning SeqRecord objects.
	"""
	for record in records:
		recFlag = True
		for phase,regexPhase in dictPhase.items():
			if regexPhase.match(str(record.seq[2:6])):
				record = record[dictTrimPhase[phase]:]
				recFlag = False
				yield record
				break
		if recFlag: 
			yield record
			continue

originalReads = SeqIO.parse(fastq1open, "fastq")
trimmedReads = trim_adaptors(originalReads,dictPhase,dictTrimPhase)
#Write reads to second command line input was gz file. 
with bgzf.BgzfWriter(sys.argv[3], "wb") as outgz:
    count1 = SeqIO.write(sequences=trimmedReads, handle=outgz, format="fastq")

with bgzf.BgzfWriter(sys.argv[4], "wb") as outgz:
    count2 = SeqIO.write(sequences=SeqIO.parse(fastq1open, "fastq"), handle=outgz, format="fastq")

print "Total sequences, fastq file 1: %s"%(count1)
print "Total sequences, fastq file 2: %s"%(count2)
