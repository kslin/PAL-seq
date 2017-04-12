import gzip
from optparse import OptionParser

import numpy as np
import pandas as pd

import helpers


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--f1", dest="FASTQ1", help="fastq file for 5' reads")
    parser.add_option("--l1", dest="LEN1", type="int", help="length of read1 used")
    parser.add_option("--f2", dest="FASTQ2", help="fastq file for 3' reads")
    parser.add_option("--l2", dest="LEN2", type="int", help="length of read2 used")
    parser.add_option("-i", "--intensity", dest="INTENSITY", help="intensity file for read2")
    parser.add_option("-o","--outfile", dest="OUTFILE", help="output file for signals")

    (options, args) = parser.parse_args()

    # read in fastq1 and make a dictionary of reads
    FASTQ1 = options.FASTQ1
    options.LEN1 = 50
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

    # read in fastq2 and make a dictionary of tail starts
    start_dict = {}
    with gzip.open(options.FASTQ2,'rb') as f:
        for line in f:
            line = line.decode()
            if line[0] == '@':
                line = line.replace('\n','').split(':')
                seq_id = ':'.join(line[2:5]).split('#')[0]
                start_dict[seq_id] = int(line[-1])


    # read in intensity file
    seq_ids = []
    t_signals = []
    starts = []
    missed1, missed2 = 0,0
    with gzip.open(options.INTENSITY,'rb') as f:
        for i, line in enumerate(f):
            line = line.decode().replace('\n','').split('\t')
            seq_id = ':'.join(line[:3])
            signal = np.array([[int(i) for i in x.split()] for x in line[4:]])
            
            # skip if you find a bunch of zeros
            if list(np.std(signal, axis=1)).count(0) > 5:
                missed1 += 1
                continue
        
            # normalize by the average intensities for each nucleotide in read1
            read1_sequence = read_dict[seq_id]
            # signal = helpers.get_normalized_intensities(signal, read1_sequence, options.LEN1, options.LEN2)
            signal = helpers.get_normalized_intensities2(signal, read1_sequence, options.LEN1, options.LEN2)

            if signal is None:
                continue

            if seq_id not in start_dict:
                continue

            # calculate t_signal from normalized intensities
            # t_signal = helpers.get_t_signal(signal)
            t_signal = helpers.get_t_signal2(signal, options.LEN1, options.LEN2)

            t_signals.append(t_signal)
            seq_ids.append(seq_id)
            starts.append(start_dict[seq_id])

    print(missed1, missed2)

    # write data as a dataframe to a file
    t_signals = pd.DataFrame(np.array(t_signals))
    t_signals['ID'] = seq_ids
    t_signals['Tail_start'] = starts
    t_signals.to_csv(options.OUTFILE,sep='\t',index=False)
