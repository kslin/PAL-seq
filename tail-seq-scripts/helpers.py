import math
import sys
import time

import numpy as np
from scipy.stats import norm

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
        intensities = intensities[config.NUM_SKIP:,:]
        read1_sequence = read1_sequence[config.NUM_SKIP:]
    except:
        print("Intensities and sequences for read1 must be longer than config.NUM_SKIP")
        sys.exit()

    assert(len(intensities) == (len1 + len2 - config.NUM_SKIP)), "Intensities length not equal to l1 and l2"
    assert(len(read1_sequence) == len1 - config.NUM_SKIP), "Read1 length must equal l1"

    # convert read1_sequence into one-hot encoding of 4 bits
    read1_sequence = np.array([[float(nt == x) for x in config.NUC_ORDER] for nt in read1_sequence])

    # assert(np.sum(read1_sequence) == (len1 - config.NUM_SKIP)), "Non-{} bases found in read1".format(config.NUC_ORDER)

    # add up the counts for each nucleotide
    read1_nt_counts = np.sum(read1_sequence, axis=0)

    # return None if one or more of the nucleotides doesn't show up in the read1_sequence
    if np.min(read1_nt_counts) == 0:
        return None

    # convert intensities to an array and split into read1 and read2 intensities
    # intensities = np.array([[int(x) for x in i.split()] for i in intensities])
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

    read2_intensities_normed = np.divide(read2_intensities, norm_vals)

    return read2_intensities_normed

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


def impute(signal, nanlimit):
    """
    Given the signal as an array, impute up to nanlimit rows of missing data.
    If there are more than nanlimit, return None
    """
    # find all the rows with all zeros
    zero_rows = (signal == [0,0,0,0]).all(axis=1)
    zero_rows = np.nonzero(zero_rows)[0]

    # if there are no rows of zeros, return signal as-is
    if len(zero_rows) == 0:
        return signal

    # if there are too many rows of zeros, return None
    elif len(zero_rows) > nanlimit:
        return None

    # otherwise, use the mean in a sliding window to fill in missing rows
    else:
        signal = signal.astype(float)
        
        for row in zero_rows:
            signal[row, :] = [np.nan]*signal.shape[1]

        new_rows = []
        for row in zero_rows:
            subrows = signal[max(0, row - nanlimit): min(len(signal), row + nanlimit), :]
            new_rows.append(np.round(np.nanmean(subrows, axis=0)))

        for i, row in enumerate(zero_rows):
            signal[row, :] = new_rows[i]

        return signal.astype(int)


def get_batch_t_signal(all_params):
    """
    Calculates normalized t-signals for a batch of sequences.
    Returns t-signal values as one long string to be written to a file.
    """
    params, len1, len2, nanlimit = all_params
    write_str = ''
    skipped = 0

    # iterate through signals in the batch
    for (seq_id, gene, start, read1, signal) in zip(*params):
        signal = signal.reshape((len1+len2, 4))

        # impute missing data
        signal = impute(signal, nanlimit)

        # skip if too much missing data
        if signal is None:
            skipped += 1
            continue

        # calculate normalization
        normed_signal = get_normalized_intensities(signal, read1, len1, len2)

        # normalize t-signal and return the output as a string
        if normed_signal is not None:
            t_signal = get_t_signal(normed_signal)
            write_str += '{}\t{}\t{}\t{}\n'.format(seq_id, gene, str(start),
                                                  '\t'.join([str(x) for x in t_signal]))

        else:
            skipped += 1
    
    return skipped, write_str


def smooth(x, window_len=21):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal. (from scipy cookbook)
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]

    return np.array([np.median(s[i-(window_len-1):i+(window_len-1)]) for i in range((window_len-1),len(x)+(window_len-1))])


def get_batch_t_signal_smooth(all_params):
    """
    Calculates normalized t-signals for a batch of sequences.
    Returns t-signal values as one long string to be written to a file.
    """
    params, len1, len2, nanlimit = all_params
    write_str = ''
    skipped = 0

    # iterate through signals in the batch
    for (seq_id, gene, start, read1, signal) in zip(*params):
        signal = signal.reshape((len1+len2, 4))

        # impute missing data
        signal = impute(signal, nanlimit)

        # skip if too much missing data
        if signal is None:
            skipped += 1
            continue

        # calculate normalization
        normed_signal = get_normalized_intensities(signal, read1, len1, len2)

        # normalize t-signal and return the output as a string
        if normed_signal is not None:
            for i in range(4):
                normed_signal[:, i] = smooth(normed_signal[:, i])

            t_signal = np.log(normed_signal[:,3])
            other_signal = np.log(np.sum(normed_signal[:,:3], axis=1))

            t_signal = t_signal - other_signal
            write_str += '{}\t{}\t{}\t{}\n'.format(seq_id, gene, str(start),
                                                  '\t'.join([str(x) for x in t_signal]))

        else:
            skipped += 1
    
    return skipped, write_str


def get_tail_length_from_emissions(emit):
    emit = list((np.array(emit) == 0) + (np.array(emit) == 1)) + [0]
    if 1 not in emit:
        return 0

    start = emit.index(1)
    end  = emit.index(0, start)
    return end - start
