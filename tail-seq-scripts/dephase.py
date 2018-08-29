from Bio import SeqIO
import sys
import tarfile

fastq1TAR=tarfile.open(name=sys.argv[1], mode='r:gz')
fastq1open=fastq1TAR.extractfile(fastq1TAR.next())
for line in fastq1open:
	print(line)
	break
for record in SeqIO.parse(fastq1open, "fastq") :
       print(record.id)
       sys.exit()
