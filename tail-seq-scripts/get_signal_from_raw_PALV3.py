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
    parser.add_option("--strand", dest="STRAND", help="Library must be stranded.") #added by TJE 20180827
    parser.add_option("-b", action="store_true", dest="ANNO_TYPE", default = False) #added by TJE 20190929 

    (options, args) = parser.parse_args()

    # if the necesary input file doesn't exist, quit
    bedfile = os.path.join(options.OUTDIR, 'read1.bed')
    if not os.path.exists(bedfile):
        raise ValueError("{} does not exist. Run intersectBed first with the same output directory.".format(bedfile))
        
    # start writing log file
    logfile = open(os.path.join(options.OUTDIR, 'logfile.txt'), 'w')

    # start timing
    t0 = time.time()

    # dedup bed input and write dropped reads
    if not options.ANNO_TYPE: bed_output = pd.read_csv(bedfile, sep='\t', header=None, usecols=[0,1,2,3,18,20], engine='c')
    if options.ANNO_TYPE: bed_output = pd.read_csv(bedfile, sep='\t', header=None, usecols=[0,1,2,3,5,15], engine='c')
    bed_output[3] = [config.fastq_header_to_ID(x) for x in bed_output[3]]
    reads_dedup, dropped_reads = preprocess_helpers.dedup_bed(bed_output, options.ANNO_TYPE)

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
    else:
        if options.STRAND == "S":
            standards = pd.read_csv(options.STANDARDS, sep='\t', header=None)
            standard_dict = {preprocess_helpers.reverse_complement(x):(y,str(z)) for (x,y,z) in zip(standards[0], standards[1],standards[2])}
        elif options.STRAND == "s": #added parsing for standard strandedness here
            standards = pd.read_csv(options.STANDARDS, sep='\t', header=None)
            standard_dict = {x:(y,str(z)) for (x,y,z) in zip(standards[0], standards[1], standards[2])}
        else:
            raise ValueError("Strandedness flag must be S or s.")

    # iterate through read2, separate standards, extract sequences
    if config.FASTQ_GZIP==True:
        fastq2open=gzip.open(options.FASTQ2, mode='rb', encoding='utf-8')
    elif config.FASTQ_GZIP==False:
        fastq2open=open(options.FASTQ2,'r')

    
    ##This softclipping dict situation needs to be updated. 
    # softClippingDict = preprocess_helpers.parse_read2_BAM(options.OUTDIR)
    # logfile.write('Length of paired end mapping reads, filtered:\t{}\n'.format(len(softClippingDict)))

    keep_dict, standard_reads = preprocess_helpers.parse_read2(fastq2open, keep_dict, standard_dict)
    fastq2open.close()
    # pdb.set_trace()

    # write table of read_IDs, 3p ends, accession IDs
    all_reads = pd.concat([standard_reads, reads_dedup], axis=0).drop_duplicates(subset=['read_ID'], keep='first')
    all_reads.to_csv(os.path.join(options.OUTDIR, 'all_read_info.txt'), sep='\t', index=False)

    logfile.write('Time to exclude reads that fail to map and extract standards\t{}\n'.format(str(datetime.timedelta(seconds=int(time.time()-t0)))))
    logfile.write('Reads found from standards\t{}\n'.format(len(standard_reads)))
    logfile.write('Reads passing mapping filter\t{}\n'.format(len(keep_dict)))
    t0 = time.time()

    # read in read2, separate short tails
    if config.FASTQ_GZIP==True:
        fastq1open=gzip.open(options.FASTQ1, mode='rb', encoding='utf-8')
    elif config.FASTQ_GZIP==False:
        fastq1open=open(options.FASTQ1,'r')

    keep_dict, standard_keep_dict, dropped_read1, num_short_tails = preprocess_helpers.parse_read1(fastq1open, keep_dict, options.OUTDIR, standard_reads, config.QUAL)

    fastq2open.close()

        
    # standard_reads[standard_reads.read_ID=='1101:9029:11191']['accession'].to_string(index = False)
    
    for read, reason in dropped_read1:
        dropped_reads_outfile.write('{}\t{}\n'.format(read, reason))

    logfile.write('Time to parse fastq2 and manually call short tails\t{}\n'.format(str(datetime.timedelta(seconds=int(time.time()-t0)))))
    logfile.write('Number short tail reads called manually\t{}\n'.format(num_short_tails))
    logfile.write('Skipped due to low quality read2:\t{}\n'.format(len(dropped_read1)))
    logfile.write('Reads for calculating normalized T-signal\t{}\n'.format(len(keep_dict)))
    t0 = time.time()
    dropped_intensity, num_reads_kept = get_signal_helpers.calculate_intensities(options.INTENSITY, keep_dict, options.OUTDIR, config.FUTURES)
    dropped_intensity_std, num_reads_kept_std = get_signal_helpers.calculate_intensities(options.INTENSITY, standard_keep_dict, options.OUTDIR, config.FUTURES, std = True)

    logfile.write('Time to calculate normalized T-signal\t{}\n'.format(str(datetime.timedelta(seconds=int(time.time()-t0)))))
    logfile.write('Skipped due to low quality intensity values\t{}\n'.format(len(dropped_intensity)))
    logfile.write('Reads for HMM\t{}\n'.format(num_reads_kept))

    for r in dropped_intensity:
        dropped_reads_outfile.write('{}\tlow_qual_intensity_vals\n'.format(r))

    # close outfiles
    dropped_reads_outfile.close()
    logfile.close()
