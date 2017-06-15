import concurrent.futures
import gzip
from optparse import OptionParser
import os
import time

import numpy as np
import pandas as pd

import config
import helpers


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--f1", dest="FASTQ1", help="fastq file for 5' reads")
    parser.add_option("--l1", dest="LEN1", type="int", help="length of read1 used")
    parser.add_option("--f2", dest="FASTQ2", help="fastq file for 3' reads")
    parser.add_option("--l2", dest="LEN2", type="int", help="length of read2 used")
    parser.add_option("-a", dest="ALIGN", help="file of alignments")
    parser.add_option("-i", "--intensity", dest="INTENSITY", help="intensity file for read2")
    parser.add_option("-o","--outdir", dest="OUTDIR", help="output directory")
    parser.add_option("-f", "--futures", dest="FUTURES", type="int", default=0, help="number of threads to use")

    (options, args) = parser.parse_args()

    # if the output directory doesn't exist, make it
    if not os.path.exists(options.OUTDIR):
        os.makedirs(options.OUTDIR)

    # read in fastq1 and make a dictionary of reads
    print("Reading fastq1...")
    FASTQ1 = options.FASTQ1
    read_dict = {}
    with gzip.open(FASTQ1, 'r') as f:
        identifier = None
        for i, line in enumerate(f):
            if (i % 4) == 0:
                identifier = helpers.get_identifier(line.decode())
            elif (i % 4) == 1:
                assert(len(line) == (options.LEN1 + 1)), "Reads in fastq1 must match given length"
                read_dict[identifier] = line.decode().replace('\n','')
            else:
                continue


    # read in alignment file and make a dictionary of reads to gene names
    print("Reading alignment file...")
    gene_dict = {}
    with open(options.ALIGN, 'rb') as f:
        for line in f:
            line = line.decode().split()
            seq_id = helpers.get_identifier(line[0])
            gene_dict[seq_id] = line[2]


    # read in fastq2 and make a dictionary of tail starts
    print("Reading fastq2...")
    start_dict = {}
    with gzip.open(options.FASTQ2,'rb') as f:
        for line in f:
            line = line.decode()
            if line[0] == '@':
                line = line.replace('\n','').split(':')
                seq_id = ':'.join(line[2:5]).split('#')[0]
                start_dict[seq_id] = int(line[-1])


    ### Read in intensity file and calculate normalized t-signal ###
    print("Writing normalized t-signal outputs to {}".format(options.OUTDIR))

    # keep track of how many reads we skip due having too many 0's
    skipped = 0
    skipped_notfound = 0

    # create an iterator that read in the data in chunks
    chunk_iterator = pd.read_csv(options.INTENSITY, delim_whitespace=True, header=None,
                                 iterator=True, chunksize=config.CHUNKSIZE)

    # these are the columns with the signal values
    signal_columns = np.arange(4*(options.LEN1 + options.LEN2)) + 4

    # if indicated, run parallel version
    if options.FUTURES > 0:
        print("Running parallel version")
        
        chunks = []

        # iterate through chunks, add gene, read, etc information, append to list of chunks
        with open(os.path.join(options.OUTDIR, 'normalized_t_signal.txt'), 'w') as outfile:
            for chunknum, chunk in enumerate(chunk_iterator):
                chunk['ID'] = chunk[0].astype(str) + ':' + chunk[1].astype(str) + ':' + chunk[2].astype(str)
                chunk['gene'] = [gene_dict[seq_id] if seq_id in gene_dict else np.nan for seq_id in chunk['ID']]
                chunk['start'] = [start_dict[seq_id] if seq_id in start_dict else np.nan for seq_id in chunk['ID']]
                chunk['read1'] = [read_dict[seq_id] if seq_id in read_dict else np.nan for seq_id in chunk['ID']]
                
                chunk = chunk.dropna()

                # record how many ids were not in the alignment/read files
                skipped_notfound += (config.CHUNKSIZE - len(chunk))

                chunks.append(((chunk['ID'].values, chunk['gene'].values, chunk['start'].values,
                                chunk['read1'].values, chunk[signal_columns].values),
                               options.LEN1, options.LEN2, config.NAN_LIMIT))

                # once we've read enough chunks, calculate t-signals in parallel
                if len(chunks) == options.FUTURES:
                    t0 = time.time()
                    full_write_str = ''
                    print("Processing up to line {}...".format((chunknum+1)*config.CHUNKSIZE))
                    with concurrent.futures.ProcessPoolExecutor() as executor:
                        results = executor.map(helpers.get_batch_t_signal, chunks)
                        # results = executor.map(helpers.get_batch_t_signal_smooth, chunks)

                        # record how many signals were skipped and write results to a file
                        for sk, wrstr in results:
                            time.sleep(0.0001)
                            skipped += sk
                            full_write_str += wrstr

                    outfile.write(full_write_str)

                    full_write_str = ''
                    chunks = []

                    print("{} seconds".format(time.time() - t0))

            # calculate t-signals for the last chunk
            with concurrent.futures.ProcessPoolExecutor() as executor:
                # results = executor.map(helpers.get_batch_t_signal, chunks)
                results = executor.map(helpers.get_batch_t_signal_smooth, chunks)
                for sk, wrstr in results:
                    skipped += sk
                    outfile.write(wrstr)

    # otherwise, run sequentially
    else:
        print("Running non-parallel version")
        
        # iterate through chunks, add gene, read, etc information
        with open(os.path.join(options.OUTDIR, 'normalized_t_signal.txt'), 'w') as outfile:
            for chunknum, chunk in enumerate(chunk_iterator):
                chunk['ID'] = chunk[0].astype(str) + ':' + chunk[1].astype(str) + ':' + chunk[2].astype(str)
                chunk['gene'] = [gene_dict[seq_id] if seq_id in gene_dict else np.nan for seq_id in chunk['ID']]
                chunk['start'] = [start_dict[seq_id] if seq_id in start_dict else np.nan for seq_id in chunk['ID']]
                chunk['read1'] = [read_dict[seq_id] if seq_id in read_dict else np.nan for seq_id in chunk['ID']]
                
                chunk = chunk.dropna()
                print(len(chunk))

                # record how many ids were not in the alignment/read files
                skipped_notfound += (config.CHUNKSIZE - len(chunk))

                # calculate t-signal
                skipped, writestr = helpers.get_batch_t_signal(((chunk['ID'].values, chunk['gene'].values, chunk['start'].values,
                                                       chunk['read1'].values, chunk[signal_columns].values),
                                                       options.LEN1, options.LEN2, config.NAN_LIMIT))

                # skipped, writestr = helpers.get_batch_t_signal_smooth(((chunk['ID'].values, chunk['gene'].values, chunk['start'].values,
                #                                        chunk['read1'].values, chunk[signal_columns].values),
                #                                        options.LEN1, options.LEN2, config.NAN_LIMIT))

                # write to file
                outfile.write(writestr)

                if ((chunknum+1) % 10) == 0:
                    print("Processed {} lines".format((chunknum+1)*config.CHUNKSIZE))


    print("Skipped due to low quality: {}".format(skipped))
    print("Skipped due to not being found in the alignment/read files: {}".format(skipped_notfound))
