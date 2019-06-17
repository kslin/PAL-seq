import datetime
from optparse import OptionParser
import os
import pickle
import time
import sys
import tarfile
import gzip
import pdb
import numpy as np
import pandas as pd

import config, preprocess_helpers, get_signal_helpers


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--f1", dest="FASTQ1", help="gzipped fastq file for 5' reads")
    parser.add_option("--f2", dest="FASTQ2", help="gzipped fastq file for 3' reads")
    parser.add_option("-i", "--intensity", dest="INTENSITY", help="gzipped intensity file for read2")
    parser.add_option("-s","--standards", dest="STANDARDS", default=None, help="file with standard sequences")
    parser.add_option("-o","--outdir", dest="OUTDIR", help="output directory")
    parser.add_option("--strand", dest="STRAND", help="strand") #added by TJE 20180827

    (options, args) = parser.parse_args()

    # if the necesary input file doesn't exist, quit
    bedfile = os.path.join(options.OUTDIR, 'read1.bed')
    if not os.path.exists(bedfile):
        print("{} does not exist. Run intersectBed first with the same output directory.".format(f))
        sys.exit()

    # start writing log file
    logfile = open(os.path.join(options.OUTDIR, 'logfile.txt'), 'w')

    # start timing
    t0 = time.time()

    # dedup bed input and write dropped reads
    bed_output = pd.read_csv(bedfile, sep='\t', header=None, usecols=[0,1,2,3,18,20], engine='c')
    bed_output[3] = [config.fastq_header_to_ID(x) for x in bed_output[3]]
    reads_dedup, dropped_reads = preprocess_helpers.dedup_bed(bed_output)

    # start writing dropped reads
    dropped_reads_outfile = open(os.path.join(options.OUTDIR, 'dropped_reads.txt'), 'w')
    dropped_reads_outfile.write('read_ID\treason\n')
    for r in dropped_reads:
        dropped_reads_outfile.write('{}\tmulti_accession_map\n'.format(r))

    # make a dictionary for reads that will continue to HMM
    keep_dict = {x: None for x in reads_dedup['read_ID']}

    logfile.write('Time to dedup BED and extract 3` ends\t{}\n'.format(str(datetime.timedelta(seconds=int(time.time()-t0)))))
    logfile.write('Skipped due to mapping to multiple accessions\t{}\n'.format(len(dropped_reads)))
    logfile.write('Reads passing dedup filter\t{}\n'.format(len(keep_dict)))
    t0 = time.time()

    # read in standards
    if options.STANDARDS is None:
        standard_dict = {}
    elif options.STRAND == "S":
        standards = pd.read_csv(options.STANDARDS, sep='\t', header=None)
        standard_dict = {preprocess_helpers.reverse_complement(x):y for (x,y) in zip(standards[0], standards[1])}
    elif options.STRAND == "s": #added parsing for standard strandedness here
        standards = pd.read_csv(options.STANDARDS, sep='\t', header=None)
        standard_dict = {x:y for (x,y) in zip(standards[0], standards[1])}
    else:
        logfile.write("No strandedness flag found.")


    # iterate through read1, separate standards, extract sequences
    if config.FASTQ_GZIP==True:
        fastq1open=gzip.open(options.FASTQ1, mode='rb')
    elif config.FASTQ_GZIP==False:
        fastq1open=open(name=options.FASTQ1,mode='r')
    else:
        fastq1Tarfile=tarfile.open(name=options.FASTQ1, mode='r:gz')
        fastq1open=fastq1Tarfile.extractfile(fastq1Tarfile.next())
    

    keep_dict, standard_reads = preprocess_helpers.parse_read1(fastq1open, keep_dict, standard_dict)
    fastq1open.close()

    # write table of read_IDs, 3p ends, accession IDs
    all_reads = pd.concat([standard_reads, reads_dedup], axis=0).drop_duplicates(subset=['read_ID'], keep='first')
    all_reads.to_csv(os.path.join(options.OUTDIR, 'all_read_info.txt'), sep='\t', index=False)

    logfile.write('Time to exclude reads that fail to map and extract standards\t{}\n'.format(str(datetime.timedelta(seconds=int(time.time()-t0)))))
    logfile.write('Reads found from standards\t{}\n'.format(len(standard_reads)))
    logfile.write('Reads passing mapping filter\t{}\n'.format(len(keep_dict)))
    t0 = time.time()

    # read in read2, separate short tails
    if config.FASTQ_GZIP==True:
        fastq2open=gzip.open(options.FASTQ2,mode='rb')
    elif config.FASTQ_GZIP==False:
        fastq2open=open(name=options.FASTQ2,mode='r')
    else:
        fastq2Tarfile=tarfile.open(name=options.FASTQ2, mode='r:gz')
        fastq2open=fastq2Tarfile.extractfile(fastq2Tarfile.next())

    softClippingDict = preprocess_helpers.parse_read2_BAM(options.OUTDIR)
    keep_dict, dropped_read2, num_short_tails = preprocess_helpers.parse_read2(fastq2open, keep_dict, options.OUTDIR, softClippingDict, standard_reads, config.QUAL)
    fastq2open.close()


    for read, reason in dropped_read2:
        dropped_reads_outfile.write('{}\t{}\n'.format(read, reason))

    logfile.write('Time to parse fastq2 and manually call short tails\t{}\n'.format(str(datetime.timedelta(seconds=int(time.time()-t0)))))
    logfile.write('Number short tail reads called manually\t{}\n'.format(num_short_tails))
    logfile.write('Skipped due to low quality read2:\t{}\n'.format(len(dropped_read2)))
    logfile.write('Reads for calculating normalized T-signal\t{}\n'.format(len(keep_dict)))
    t0 = time.time()
    

    # keep_dict = {'1101:9021:2000#TAGTGC': (preprocess_helpers.reverse_complement('ACCAAAAATCTGTCACAGAATTTTGAGACCATTAAAACAAGTTTAATGAN'), 0),
    #              '1101:9793:1998#TAGTGC': (preprocess_helpers.reverse_complement('GGCTGGCCTGTACACTGACTTGAGACCAATAAAAGTGCACACCTTACCTN'), 0),
    #              '1101:16455:1995#TAGTGC': (preprocess_helpers.reverse_complement('CCCTAAAATTGGTTTCAAGCCAATCTCATATCCTATATGTCTTTCTCAAN'), 0), 
    #              '1101:16631:1996#TAGTGC': (preprocess_helpers.reverse_complement('CGGCTGTGGGAATGAATCATTGAAGTAATAAACTACAGTGGTTGATCCAN'), 0)}


    # pickle.dump(keep_dict, open(os.path.join(options.OUTDIR, 'keep_dict.pickle'), "w"))
    # keep_dict = pickle.load(open(os.path.join(options.OUTDIR, 'keep_dict.pickle'), "rb"))   

    dropped_intensity, num_reads_kept = get_signal_helpers.calculate_intensities(options.INTENSITY, keep_dict, options.OUTDIR, config.FUTURES)
    
    logfile.write('Time to calculate normalized T-signal\t{}\n'.format(str(datetime.timedelta(seconds=int(time.time()-t0)))))
    logfile.write('Skipped due to low quality intensity values\t{}\n'.format(len(dropped_intensity)))
    logfile.write('Reads for HMM\t{}\n'.format(num_reads_kept))

    for r in dropped_intensity:
        dropped_reads_outfile.write('{}\tlow_qual_intensity_vals\n'.format(r))

    # close outfiles
    dropped_reads_outfile.close()
    logfile.close()
