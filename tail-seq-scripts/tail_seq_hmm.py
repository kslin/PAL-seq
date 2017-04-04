import gzip
from optparse import OptionParser
import sys

import pandas as pd

import helpers


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--f1", dest="FASTQ1", help="fastq file for 5' reads")
    parser.add_option("--l1", dest="LEN1", type="int", help="length of read1 used")
    parser.add_option("--f2", dest="FASTQ2", help="fastq file for 3' reads")
    parser.add_option("--l2", dest="LEN2", type="int", help="length of read2 used")
    parser.add_option("-m", "--mapping", dest="MAPPING", help="mapping file for fastq1")
    parser.add_option("-i", "--intensity", dest="INTENSITY", help="intensity file for read2")
    parser.add_option("-t","--trainsize", default=10000, type="int", dest="TRAINSIZE", help="length of training set")
    parser.add_option("-o","--outdir", dest="OUTDIR", help="directory of outputs")
    parser.add_option("-s", "--standards", dest="STANDARDS", default=None,
                      help="file with names of standards")

    (options, args) = parser.parse_args()

    # read in fastq1 and make a dictionary of reads
    fastq1_dict = {}
    with gzip.open(options.FASTQ1, 'rb') as f:
        identifier = None
        for i, line in enumerate(f):
            if (i % 4) == 0:
                identifier = helpers.get_identifier(line.decode())
            elif (i % 4) == 1:
                assert(len(line) == (options.LEN1 + 1)), 'Reads in fastq1 must match given length'
                fastq1_dict[identifier] = line.decode().replace('\n','')
            else:
                continue

    print(len(fastq1_dict))