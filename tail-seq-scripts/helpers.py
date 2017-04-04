import math
import sys

import matplotlib.pyplot as plt
import numpy as np

import config


def get_identifier(raw_header):
    return ':'.join(raw_header.split(':')[2:5]).split('#')[0]

def fix_vals(x):
    if x > 0:
        return x
    return 1.0

v_fix_vals = np.vectorize(fix_vals)

def get_normalized_intensities(intensities, read1_sequence, len1, len2):
    """Get average intensities for each nucleotide for use in normalizing"""

    # skip the first 4 elements, which identify the read
    # then skip the number of starting nucleotides specified in the config file
    try:
        intensities = intensities.replace('\n','').split('\t')[(4 + config.NUM_SKIP):]
        read1_sequence = read1_sequence[config.NUM_SKIP:]
    except:
        print("Intensities and sequences for read1 must be longer than config.NUM_SKIP")
        sys.exit()

    assert(len(intensities) == (len1 + len2 - config.NUM_SKIP)), "Intensities length not equal to l1 and l2"
    assert(len(read1_sequence) == len1 - config.NUM_SKIP), "Read1 length must equal l1"

    # convert read1_sequence into one-hot encoding of 4 bits
    read1_sequence = np.array([[int(nt == x) for x in config.NUC_ORDER] for nt in read1_sequence])

    assert(np.sum(read1_sequence) == (len1 - config.NUM_SKIP)), "Non-(A,C,G,T) bases found in read1"

    # add up the counts for each nucleotide
    read1_nt_counts = np.sum(read1_sequence, axis=0)

    # return None if one or more of the nucleotides doesn't show up in the read1_sequence
    if np.min(read1_nt_counts) == 0:
        print("Warning: ")
        return None

    # convert intensities to an array and split into read1 and read2 intensities
    intensities = np.array([[int(x) for x in i.split()] for i in intensities])
    read1_intensities = intensities[:len1 - config.NUM_SKIP]
    read2_intensities = intensities[-1 * len2:]

    # if any intensities are negative, convert to 1
    read1_intensities = v_fix_vals(read1_intensities)
    read2_intensities = v_fix_vals(read2_intensities)

    # multiply read1 intensities by one-hot to get intensities for the base called
    read1_intensities = np.multiply(read1_sequence, read1_intensities)

    # add up read1_intensities for each nucleotide
    read1_intensities = np.sum(read1_intensities, axis=0)

    # get average read1_intensities for each nucleotide
    norm_vals = np.divide(read1_intensities, read1_nt_counts)

    return np.divide(read2_intensities, norm_vals)


def get_t_signal(intensities, upperbound=config.UPPERBOUND, lowerbound=config.LOWERBOUND):

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
    t_signal = np.minimum(t_signal, upperbound)
    t_signal = np.maximum(t_signal, lowerbound)

    return t_signal

