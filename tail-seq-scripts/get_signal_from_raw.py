from optparse import OptionParser
import os
import time
import sys

import numpy as np
import pandas as pd

import config, preprocess_helpers, get_signal_helpers


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--f1", dest="FASTQ1", help="gzipped fastq file for 5' reads")
    parser.add_option("--f2", dest="FASTQ2", help="gzipped fastq file for 3' reads")
    parser.add_option("-b", dest="BEDFILE", help="output of intersect bed")
    parser.add_option("-i", "--intensity", dest="INTENSITY", help="gzipped intensity file for read2")
    parser.add_option("-s","--standards", dest="STANDARDS", default=None, help="file with standard sequences")
    parser.add_option("-o","--outdir", dest="OUTDIR", help="output directory")
    parser.add_option("-t","--short_tails", dest="SHORT_TAILS", help="output file for manually called short tails")
    parser.add_option("-d","--dropped", dest="DROPPED", help="output file for dropped reads")
    parser.add_option("-f", "--futures", dest="FUTURES", type="int", default=1, help="number of threads to use")

    (options, args) = parser.parse_args()

    # if the output directory doesn't exist, make it
    if not os.path.exists(options.OUTDIR):
        os.makedirs(options.OUTDIR)

    ## dedup bed input and write dropped reads
    bedfile = pd.read_csv(options.BEDFILE, sep='\t', header=None, usecols=[0,1,2,3,18,20])
    bedfile[3] = [config.fastq_header_to_ID(x) for x in bedfile[3]]
    reads_dedup, dropped_reads = preprocess_helpers.dedup_bed(bedfile)

    # start writing dropped reads
    dropped_reads_outfile = open(options.DROPPED, 'w')
    dropped_reads_outfile.write('read_ID\treason\n')
    for r in dropped_reads:
        dropped_reads_outfile.write('{}\tmulti_gene_map\n'.format(r))

    # make a dictionary for reads that will continue to HMM
    keep_dict = {x: None for x in reads_dedup['read_ID']}

    print('Skipped due to mapping to multiple genes: {}'.format(len(dropped_reads)))
    print('Reads moving forward: {}'.format(len(keep_dict)))

    # read in standards
    if options.STANDARDS is None:
        standard_dict = {}
    else:
        standards = pd.read_csv(options.STANDARDS, sep='\t', header=None)
        standard_dict = {x:y for (x,y) in zip(standards[0], standards[1])}

    # iterate through read1, separate standards, extract sequences
    keep_dict, standard_reads = preprocess_helpers.parse_read1(options.FASTQ1, keep_dict, standard_dict)
    print('Found {} reads from standards'.format(len(standard_reads)))
    print('Reads moving forward: {}'.format(len(keep_dict)))

    all_reads = pd.concat([standard_reads, reads_dedup], axis=0).drop_duplicates(subset=['read_ID'], keep='first')
    all_reads.to_csv(os.path.join(options.OUTDIR, 'all_read_info.txt'), sep='\t', index=False)

    print(len(keep_dict))

    # read in read2, separate short tails
    keep_dict, dropped_read2 = preprocess_helpers.parse_read2(options.FASTQ2, keep_dict, options.SHORT_TAILS)

    print('Skipped due to low quality read2: {}'.format(len(dropped_read2)))
    print('Reads moving forward: {}'.format(len(keep_dict)))
    for read, reason in dropped_read2:
        dropped_reads_outfile.write('{}\t{}\n'.format(read, reason))

    dropped_intensity, num_reads_kept = get_signal_helpers.calculate_intensities(options.INTENSITY, keep_dict, options.OUTDIR, options.FUTURES)

    print("Skipped due to low quality intensity values: {}".format(len(dropped_intensity)))
    print('Reads moving forward: {}'.format(num_reads_kept))
    for r in dropped_intensity:
        dropped_reads_outfile.write('{}\tlow_qual_intensity_vals\n'.format(r))

    # close dropped reads file
    dropped_reads_outfile.close()
