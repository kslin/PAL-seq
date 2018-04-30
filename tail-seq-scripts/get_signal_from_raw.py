from optparse import OptionParser
import os
import pickle
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

    print('Deduplicating BED file and extracting 3` ends...')
    t0 = time.time()

    # dedup bed input and write dropped reads
    bedfile = pd.read_csv(options.BEDFILE, sep='\t', header=None, usecols=[0,1,2,3,18,20], engine='c')
    bedfile[3] = [config.fastq_header_to_ID(x) for x in bedfile[3]]
    reads_dedup, dropped_reads = preprocess_helpers.dedup_bed(bedfile)

    # start writing dropped reads
    dropped_reads_outfile = open(options.DROPPED, 'w')
    dropped_reads_outfile.write('read_ID\treason\n')
    for r in dropped_reads:
        dropped_reads_outfile.write('{}\tmulti_gene_map\n'.format(r))

    # make a dictionary for reads that will continue to HMM
    keep_dict = {x: None for x in reads_dedup['read_ID']}

    print('{:.3f} seconds'.format(time.time() - t0))
    print('Skipped due to mapping to multiple genes: {}'.format(len(dropped_reads)))
    print('Reads moving forward: {}'.format(len(keep_dict)))

    print('Excluding reads that fail to map and extracting standards...')
    t0 = time.time()

    # read in standards
    if options.STANDARDS is None:
        standard_dict = {}
    else:
        standards = pd.read_csv(options.STANDARDS, sep='\t', header=None)
        standard_dict = {preprocess_helpers.reverse_complement(x):y for (x,y) in zip(standards[0], standards[1])}

    # iterate through read1, separate standards, extract sequences
    keep_dict, standard_reads = preprocess_helpers.parse_read1(options.FASTQ1, keep_dict, standard_dict)

    # write table of read_IDs, 3p ends, gene IDs
    all_reads = pd.concat([standard_reads, reads_dedup], axis=0).drop_duplicates(subset=['read_ID'], keep='first')
    all_reads.to_csv(os.path.join(options.OUTDIR, 'all_read_info.txt'), sep='\t', index=False)

    print('{:.3f} seconds'.format(time.time() - t0))
    print('Found {} reads from standards'.format(len(standard_reads)))
    print('Reads moving forward: {}'.format(len(keep_dict)))

    print('Calling short tails manually...')
    t0 = time.time()

    # read in read2, separate short tails
    keep_dict, dropped_read2 = preprocess_helpers.parse_read2(options.FASTQ2, keep_dict, options.SHORT_TAILS)

    for read, reason in dropped_read2:
        dropped_reads_outfile.write('{}\t{}\n'.format(read, reason))

    print('{:.3f} seconds'.format(time.time() - t0))
    print('Skipped due to low quality read2: {}'.format(len(dropped_read2)))
    print('Reads moving forward: {}'.format(len(keep_dict)))

    print('Calculating normalized T-signal...')
    t0 = time.time()
    

    # keep_dict = {'1101:9021:2000#TAGTGC': (preprocess_helpers.reverse_complement('ACCAAAAATCTGTCACAGAATTTTGAGACCATTAAAACAAGTTTAATGAN'), 0),
    #              '1101:9793:1998#TAGTGC': (preprocess_helpers.reverse_complement('GGCTGGCCTGTACACTGACTTGAGACCAATAAAAGTGCACACCTTACCTN'), 0),
    #              '1101:16455:1995#TAGTGC': (preprocess_helpers.reverse_complement('CCCTAAAATTGGTTTCAAGCCAATCTCATATCCTATATGTCTTTCTCAAN'), 0), 
    #              '1101:16631:1996#TAGTGC': (preprocess_helpers.reverse_complement('CGGCTGTGGGAATGAATCATTGAAGTAATAAACTACAGTGGTTGATCCAN'), 0)}


    # pickle.dump(keep_dict, open(os.path.join(options.OUTDIR, 'keep_dict.pickle'), "w"))
    # keep_dict = pickle.load(open(os.path.join(options.OUTDIR, 'keep_dict.pickle'), "r"))

    # new_keep_dict = {}
    # for key, val in keep_dict.items():
    #     key1,key2,key3 = key.split('#')[0].split(':')
    #     if key1 not in new_keep_dict:
    #         new_keep_dict[key1] = {key2:{key3:val}}
    #     else:
    #         if key2 not in new_keep_dict[key1]:
    #             new_keep_dict[key1][key2] = {key3:val}
    #         else:
    #             new_keep_dict[key1][key2][key3] = val

    # print(len(new_keep_dict))

    # pickle.dump(new_keep_dict, open(os.path.join(options.OUTDIR, 'new_keep_dict.pickle'), "w"))
    # new_keep_dict = pickle.load(open(os.path.join(options.OUTDIR, 'new_keep_dict.pickle'), "r"))

    dropped_intensity, num_reads_kept = get_signal_helpers.calculate_intensities(options.INTENSITY, keep_dict, options.OUTDIR, options.FUTURES)

    print('{:.3f} seconds'.format(time.time() - t0))
    print("Skipped due to low quality intensity values: {}".format(len(dropped_intensity)))
    print('Reads moving forward: {}'.format(num_reads_kept))
    for r in dropped_intensity:
        dropped_reads_outfile.write('{}\tlow_qual_intensity_vals\n'.format(r))

    # close dropped reads file
    dropped_reads_outfile.close()
