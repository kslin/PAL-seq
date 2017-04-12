import gzip
from optparse import OptionParser
import sys

from hmmlearn.hmm import GaussianHMM
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config
import helpers
import models


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--f1", dest="FASTQ1", help="fastq file for 5' reads")
    parser.add_option("--l1", dest="LEN1", type="int", help="length of read1 used")
    parser.add_option("--f2", dest="FASTQ2", help="fastq file for 3' reads")
    parser.add_option("--l2", dest="LEN2", type="int", help="length of read2 used")
    parser.add_option("-m", "--mapping", dest="MAPPING", help="mapping file for fastq1")
    parser.add_option("-i", "--intensity", dest="INTENSITY", help="intensity file for read2")
    parser.add_option("-t","--trainsize", default=1000, type="int", dest="TRAINSIZE", help="length of training set")
    parser.add_option("-o","--outdir", dest="OUTDIR", help="directory of outputs")
    parser.add_option("-s", "--standards", dest="STANDARDS", default=None,
                      help="file with names of standards")

    (options, args) = parser.parse_args()

    print("Creating read dictionary")
    read_dict = {}

    # read in fastq1 and make a dictionary of reads
    with gzip.open(options.FASTQ1, 'rb') as f:
        identifier = None
        for i, line in enumerate(f):
            if (i % 4) == 0:
                identifier = helpers.get_identifier(line.decode())
            elif (i % 4) == 1:
                assert(len(line) == (options.LEN1 + 1)), "Reads in fastq1 must match given length"
                read_dict[identifier] = [line.decode().replace('\n','')]
            else:
                continue

    # read in fastq2 and add on the position where the tail starts
    with gzip.open(options.FASTQ2, 'rb') as f:
        identifier = None
        for i, line in enumerate(f):
            if (i % 4) == 0:
                line = line.decode().replace('\n','')
                identifier = helpers.get_identifier(line)
                tail_ix = line.split('::')[1]
                if identifier not in read_dict:
                    print("Warning: {} found in fastq2 but not in fastq1".format(identifier))
                else:
                    read_dict[identifier].append(int(tail_ix))
            else:
                continue

    assert(len(read_dict) >= options.TRAINSIZE), "Data must be at least as long as the training size"

    print("Reading intensity data")
    adata, cdata, gdata, tdata = [], [], [], []
    ids = []

    # read training set from intensity file, store normalized signal for T channel
    with gzip.open(options.INTENSITY, 'rb') as f:
        for i, line in enumerate(f):
            line = line.decode().replace('\n','').split('\t')

            # extract ID 
            seq_id = ':'.join(line[:3])

            # check that ID is in read_dict
            if seq_id not in read_dict:
                print("Warning: sequence ID {} from intensity file not found in fastq1".format(seq_id))
                continue

            # get sequence and intensities
            seq_info = read_dict[seq_id]
            if len(seq_info) != 2:
                print("Warning: {} found in only one of fastq1 and fastq2".format(seq_id))
                continue

            read1_seq, tail_ix = seq_info
            intensities = line[4:]

            # get normalized intensity values for read2
            r2_normed_intensities = helpers.get_normalized_intensities(intensities, read1_seq, options.LEN1, options.LEN2)
            if r2_normed_intensities is not None:
                adata.append(r2_normed_intensities[:, 0])
                cdata.append(r2_normed_intensities[:, 1])
                gdata.append(r2_normed_intensities[:, 2])
                tdata.append(r2_normed_intensities[:, 3])
                ids.append(seq_id)

            if i == 3000:
                break

    adata = np.array(adata)
    cdata = np.array(cdata)
    gdata = np.array(gdata)
    tdata = np.array(tdata)
    ids = np.array(ids)

    # separate out training set
    train_ix = np.random.choice(range(len(adata)), size=options.TRAINSIZE)
    training_ids = ids[train_ix]

    training = adata[train_ix, :]
    training_df = pd.DataFrame(training)
    training_df['ID'] = training_ids
    training_df.to_csv('/lab/bartel4_ata/kathyl/Tail_Seq/data/atraining.csv', index=False)

    training = cdata[train_ix, :]
    training_df = pd.DataFrame(training)
    training_df['ID'] = training_ids
    training_df.to_csv('/lab/bartel4_ata/kathyl/Tail_Seq/data/ctraining.csv', index=False)

    training = gdata[train_ix, :]
    training_df = pd.DataFrame(training)
    training_df['ID'] = training_ids
    training_df.to_csv('/lab/bartel4_ata/kathyl/Tail_Seq/data/gtraining.csv', index=False)

    training = tdata[train_ix, :]
    training_df = pd.DataFrame(training)
    training_df['ID'] = training_ids
    training_df.to_csv('/lab/bartel4_ata/kathyl/Tail_Seq/data/ttraining.csv', index=False)

    # write to a file
    # np.savetxt('/lab/bartel4_ata/kathyl/Tail_Seq/data/training.csv', training, delimiter=',')

    # print(training.shape)

                    
    # Build HMM
    # model = models.Transition()











