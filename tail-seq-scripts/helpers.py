import numpy as np

import config


def get_identifier(raw_header):
    """Extract sequencing cluster identifier from header. """
    return ':'.join(raw_header.split(':')[2:5]).split('#')[0]


def fix_vals(x):
    """Convert negative values to 1."""
    if x > 0:
        return x
    return 1.0


# vectorize fix_vals
v_fix_vals = np.vectorize(fix_vals)


def get_normalized_intensities(intensities, read1_sequence, len1, len2):
    """Use average intensities for the nucleotides in read1 to normalize the signals for read2.

    Arguments:
        intensities - numpy array of shape (len1 + len2) x (number of nucleotides i.e. 4)
        read1_sequence - string sequence of read1
        len1 - length of read1 as an int
        len2 - length of read2 as an int
    """
    # skip the number of starting nucleotides specified in the config file
    try:
        intensities = intensities[config.NUM_SKIP:,:]
        read1_sequence = read1_sequence[config.NUM_SKIP:]
    except:
        print("Intensities and sequences for read1 must be longer than config.NUM_SKIP")

    if len(intensities) != (len1 + len2 - config.NUM_SKIP):
        raise ValueError("Intensities length not equal to len1 + len2")

    if len(read1_sequence) != (len1 - config.NUM_SKIP):
        raise ValueError("Read1 length must equal len1")

    # convert read1_sequence into one-hot encoding of 4 bits
    read1_sequence = np.array([[float(nt == x) for x in config.NUC_ORDER] for nt in read1_sequence])

    # add up the counts for each nucleotide
    read1_nt_counts = np.sum(read1_sequence, axis=0)

    # return None if one or more of the nucleotides doesn't show up in the read1_sequence
    if np.min(read1_nt_counts) == 0:
        return None

    # convert intensities to an array and split into read1 and read2 intensities
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

    # divide read2 intensities by normalization values
    read2_intensities_normed = np.divide(read2_intensities, norm_vals)

    return read2_intensities_normed

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
    t_signal = np.minimum(t_signal, config.UPPERBOUND)
    t_signal = np.maximum(t_signal, config.LOWERBOUND)

    return t_signal


def impute(signal, nanlimit):
    """
    Given the signal as an array, impute up to nanlimit rows of missing data.
    If there are more than nanlimit, return None.
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
    Calculate normalized t-signals for a batch of sequences.
    Return t-signal values as one long string to be written to a file.
    """
    params, len1, len2 = all_params
    write_str = ''
    skipped = 0

    # iterate through signals in the batch
    for (seq_id, gene, start, read1, signal) in zip(*params):
        signal = signal.reshape((len1+len2, 4))

        # impute missing data
        signal = impute(signal, config.NAN_LIMIT)

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


def get_tail_length_from_emissions(emit):
    """
    Get the tail length from HMM emissions.
    Return the longest string of contiguous 0's before the first 1
    """
    emit = list((np.array(emit) == 0) + (np.array(emit) == 1)) + [0]
    if 1 not in emit:
        return 0

    start = emit.index(1)
    end  = emit.index(0, start)
    return end - start
