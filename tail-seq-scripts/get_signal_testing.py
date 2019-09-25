import concurrent.futures
import io
import gzip
import os
import time

import numpy as np
import pandas as pd

import config

def get_normalized_intensities(intensities, concatSequence):
    """Use average intensities for the nucleotides in read1 to normalize the signals for read2.

    Arguments:
        intensities: numpy array of shape (config.LEN1 + config.LEN2) x (number of nucleotides i.e. 4)
        concatSequence: string sequence of read1
    """

    # skip the number of starting nucleotides specified in the config file
    try:
        # intensities = np.vstack((intensities[config.NUM_SKIP:(config.LEN1 - config.NUM_SKIP_2),:],intensities[-1 * config.LEN2:,:]))
        # concatSequence = concatSequence[config.NUM_SKIP:(config.LEN1 - config.NUM_SKIP_2)]
        intensitiesTrim = intensities[config.NUM_SKIP:(config.LEN1 + config.LEN2 - config.NUM_SKIP_2),:]
        concatSequenceTrim = concatSequence[config.NUM_SKIP:(config.LEN1 + config.LEN2 - config.NUM_SKIP_2)]
    except:
        raise ValueError("Intensities and sequences for read1 must be longer than config.NUM_SKIP")

    # if len(intensities) != (config.LEN1 + config.LEN2 - config.NUM_SKIP - config.NUM_SKIP_2):
    if len(intensitiesTrim) != (config.LEN1 + config.LEN2 - config.NUM_SKIP - config.NUM_SKIP_2):
        raise ValueError("Intensities length not equal to config.LEN1 + config.LEN2")

    # if len(concatSequence) != (config.LEN1 - config.NUM_SKIP - config.NUM_SKIP_2):
    if len(concatSequenceTrim) != (config.LEN1 + config.LEN2 - config.NUM_SKIP - config.NUM_SKIP_2):
        raise ValueError("Read1 length must equal to intensity length")

    # convert concatSequence into one-hot encoding of 4 bits
    ## Changed on 2019 06 13.

    ### Working on this line on 2019 09 23.
    concatSequenceTrim = np.array([[float(nt == x) for x in config.NUC_ORDER] for nt in concatSequenceTrim])

    # add up the counts for each nucleotide
    concatSequenceNtCounts = np.sum(concatSequenceTrim, axis=0)

    # return None if one or more of the nucleotides doesn't show up in the concatSequence
    if np.min(concatSequenceNtCounts) == 0:
        return None

    # convert intensities to an array and split into read1 and read2 intensities
    concatSequenceIntensityTrim = intensitiesTrim[:] ##I think it's this line.
    read2_intensities = intensities[-1 * config.NUM_SKIP_2:]

    # multiply read1 intensities by one-hot to get intensities for the base called
    concatSequenceIntensityTrim = np.multiply(concatSequenceTrim, concatSequenceIntensityTrim)
    # add up concatSequenceIntensityTrim for each nucleotide
    concatSequenceIntensityTrim = np.sum(concatSequenceIntensityTrim, axis=0)

    # get average concatSequenceIntensityTrim for each nucleotide
    norm_vals = np.divide(concatSequenceIntensityTrim, concatSequenceNtCounts)

    # divide read2 intensities by normalization values
    read2_intensities_normed = np.divide(read2_intensities, norm_vals)

    #Subtract the value of the background
    read2_intensities_normed = np.subtract(read2_intensities_normed[:,], read2_intensities_normed[0,])

def get_t_signal(intensities):
    """Divide T intensity by the sum of the other intensities to get the t-signal"""

    # find which index T is at in the intensity file
    t_index = config.NUC_ORDER.index('T')
    other_index = list(range(len(config.NUC_ORDER)))
    other_index.remove(t_index)
    
    # separate T channel from the others
    t_channel = intensities[:, t_index]
    other_channels = intensities[:, other_index]
    background = np.sum(other_channels, axis=1)

    t_signal = np.log2(np.divide(t_channel, background))

    # bound the signal
    # t_signal = np.minimum(t_signal, config.UPPERBOUND)
    # t_signal = np.maximum(t_signal, config.LOWERBOUND)

    return t_signal

def get_batch_t_signal(params):
    """
    Calculate normalized t-signals for a batch of sequences.
    Return t-signal values as one long string to be written to a file.
    """

    write_str = ''
    skipped = []
    # iterate through signals in the batch
    for (read_ID, read2, start, signal) in zip(*params):

        signal = signal.reshape((config.LEN1+config.LEN2, 4))
        # # impute missing data
        # signal = impute(signal, config.NAN_LIMIT)

        # # skip if too much missing data
        # if signal is None:
        #     skipped.append(read_ID)
        #     continue

        # calculate normalization
        normed_signal = get_normalized_intensities(signal, read2)
        # normalize t-signal and return the output as a string
        if normed_signal is not None:
            t_signal = get_t_signal(normed_signal)
            write_str += '{}\t{}\t{}\n'.format(read_ID, str(start),
                                                  '\t'.join(['{:.3f}'.format(x) for x in t_signal]))

        else:
            skipped.append(read_ID)
    
    return skipped, write_str

def main():
    infile = open('sample1_hits_HEAD10000.txt', 'r')
    fastq1 = open('sample1_2_HEAD40000.txt','r')
    fastq2 = open('sample1_2_HEAD40000.txt','r')

    next(fastq1)
    next(fastq2)

    fq1Seq = next(fastq1).strip()
    fq2Seq = next(fastq2).strip()

    keep_dict = {'1101:1888:1994': (fq1Seq + fq2Seq, 10)}
    line_num = 0
    ix = 0
    IDs, read2s, starts, intensity_values = [], [], [], []
    chunk = (IDs, read2s, starts, np.array(intensity_values, dtype=int))
    for line in infile:
            line = line.split()
            line_num += 1
            read_ID = '1101:1888:1994'
    
            try:
                print(read_ID)
                seq, TailBeginLength = keep_dict[read_ID]
            except:
                break
    
            IDs.append(read_ID)
            read2s.append(seq)
            starts.append(TailBeginLength)
            intensity_values.append(line[config.SIGNAL_COL_START: config.SIGNAL_COL_END])
            ix += 1
            chunk = (IDs, read2s, starts, np.array(intensity_values, dtype=int))
            # if indicated, run parallel version
            sk, write_str = get_batch_t_signal(chunk)
            break

if __name__ == "__main__":
    main()


