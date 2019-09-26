import gzip
import numpy as np
import config
import scipy.signal as ss

def read_training_set(infile, file_length, random_seed=0):
    """
    Get training data for linear model.
    Return the sequences concatenated together and a list of sequence lengths.
    """

    training_values = []
    training_seq_lengths = []
    with open(infile, 'r') as f:
        for line in f:
            if line_counter == current_ix:
                line = line[:-1].split('\t') # remove newline character and split by tab
                tail_start = int(line[1])
                line = smooth(line)
                if (ix_counter % 2) == 0:
                    training_values.append([config.START_SIGNAL + 1])
                else:
                    training_values.append([config.START_SIGNAL - 1])
                for val in line[2 + tail_start:]:
                    training_values.append([float(val)])
                training_seq_lengths.append(config.LEN2 - tail_start + 1)

                # update counters
                ix_counter += 1
                if ix_counter == train_size:
                    break

                current_ix = keep_ix[ix_counter]

            line_counter += 1

    return np.array(training_values), training_seq_lengths


def get_tail_length_from_emissions(emit, twostate):
    """
    Get the tail length from HMM emissions.
    Return the longest string of contiguous 0's before the first 1
    """
    if twostate:
        emit = list(emit == 0) + [0]
    else:
        emit = list((emit == 0) + (emit == 1)) + [0]
    if 1 not in emit:
        return 0

    start = emit.index(1)
    end  = emit.index(0, start)
    return end - start

def smooth(signal):
    """
    Modified to smooth only after the tail starts.
    TJE 2019 06 09 
    """
    tail_starts = int(signal[1])
    f64Signal = np.array(signal[(2 + tail_starts):]).astype(np.float)
    f64SignalFilt = ss.medfilt(f64Signal,11)
    strSignalFilt = np.insert(f64SignalFilt.astype(str),0,signal[0:(2 + tail_starts)])
    return strSignalFilt
