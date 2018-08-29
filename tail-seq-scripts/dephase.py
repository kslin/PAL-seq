#Opens the fastq files from their original positions on solexa public as tarfiles, then uses Biopython to parse
#Must be run with python 2.7
#python2.7 tail-seq-scripts/dephase.py /lab/solexa_public/Bartel/160311_WIGTC-HISEQA_C73V2ACXX/QualityScore/CTAGTC-s_4_1_sequence.txt.tar.gz

from Bio import SeqIO
import sys
import tarfile
import re

fastq1TAR=tarfile.open(name=sys.argv[1], mode='r:gz')
fastq1open=fastq1TAR.extractfile(fastq1TAR.next())
counter = 0

dictPhase = {'objPhase1':re.compile("A..."),'objPhase2':re.compile("GA.."),'objPhase3':re.compile("CGA."),'objPhase4':re.compile("TCGA")}
for record in SeqIO.parse(fastq1open, "fastq") :
	counter+=1
	print(record.seq[2:6])
	for phase,regexPhase in dictPhase.items():
		if regexPhase.match(str(record.seq[2:6])): print phase
	if counter == 10: sys.exit()
