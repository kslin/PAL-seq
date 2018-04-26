from optparse import OptionParser

import numpy as np
import pandas as pd


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-i", dest="INFILE", help="path to read1 fastq file")
    parser.add_option("-s","--standards", dest="STANDARDS", help="file with standard sequences")
    parser.add_option("-o","--outfile", dest="OUTFILE", help="output file")
    parser.add_option("-p","--outfile_standards", dest="OUTFILE_STANDARDS", help="output file for standards")

    (options, args) = parser.parse_args()

    # read in standards
    standards = pd.read_csv(options.STANDARDS, sep='\t', header=None)
    standard_dict = {x:y for (x,y) in zip(standards[0], standards[1])}

    # read in fastq file and pull out standards
    with open(options.OUTFILE, 'w') as outfile:
        with open(options.OUTFILE_STANDARDS, 'w') as standard_outfile:
            with open(options.INFILE, 'r') as infile:
                while True:
                    line1 = infile.readline()
                    if len(line1) == 0:
                        break

                    if line1[0] != '@':
                        raise ValueError('fastq file entries must start with @')

                    seq = infile.readline()
                    line3 = infile.readline()
                    line4 = infile.readline()

                    for standard, name in standard_dict.items():
                        if standard in seq:
                            read_ID = line1[1:].replace('\n','')
                            standard_outfile.write('\t'.join(['standard','0','0','+',standard_dict[standard], read_ID]) + '\n')
                            break

                    else:
                        outfile.write(line1 + seq + line3 + line4)

