import concurrent.futures
import io
import gzip
import os
import time

import numpy as np
import pandas as pd

import config

'''
This code contains helper functions for parsing the instensity files from an 
Illumina Hiseq 2500. For a new user, the most important thing to note is 
line 74. This line computes the signal that is used as the normalization for 
the T channel when computing the T signal. The normalization is important because 
different clusters will have different sizes and thus intensities. 

There are two at least two ways to think about this normalization:
1) Derive a constant that results from the average intensity of each base in 
read 1 using the signal for that base. E.g. if the called base in A, the signal
from the A channel is used to compute the normalization constant. 
2) Derive a constant that results from the average intensity of each base in 
read 1 using the signal for bases other than the called base. This assesses the 
background instead of the signal: i.e. if the called base is A, the signal from
C, G, and T channels are used to compute the normalization constant. 

We used normalization strategy 2 for the splint ligation runs and normalization 
strategy 1 for the direct ligation runs (toggled by uncommenting line 73). 
In practice, either might be appropriate depending on the run, but it's 
important to make that assessment depending on how the standards are called 
after the pipeline is run. 
'''

def fix_vals(x):
    """Convert negative values to 1."""
    if x > 0:
        return x
    return 1.0


# vectorize fix_vals
v_fix_vals = np.vectorize(fix_vals)


def get_normalized_intensities(intensities, read1_sequence):
    """Use average intensities for the nucleotides in read1 to normalize the signals for read2.

    Arguments:
        intensities: numpy array of shape (config.LEN1 + config.LEN2) x (number of nucleotides i.e. 4)
        read1_sequence: string sequence of read1
    """

    # skip the number of starting nucleotides specified in the config file
    try:
        intensities = np.vstack((intensities[config.NUM_SKIP:(config.LEN1 - config.NUM_SKIP_2),:],intensities[-1 * config.LEN2:,:]))
        read1_sequence = read1_sequence[config.NUM_SKIP:(config.LEN1 - config.NUM_SKIP_2)]

        # intensities = intensities[config.NUM_SKIP:,:]
        # read1_sequence = read1_sequence[config.NUM_SKIP:]
    except:
        raise ValueError("Intensities and sequences for read1 must be longer than config.NUM_SKIP")

    if len(intensities) != (config.LEN1 + config.LEN2 - config.NUM_SKIP - config.NUM_SKIP_2):
    # if len(intensities) != (config.LEN1 + config.LEN2 - config.NUM_SKIP):
        raise ValueError("Intensities length not equal to config.LEN1 + config.LEN2")

    if len(read1_sequence) != (config.LEN1 - config.NUM_SKIP - config.NUM_SKIP_2):
    # if len(read1_sequence) != (config.LEN1 - config.NUM_SKIP):
        raise ValueError("Read1 length must equal config.LEN1")

    # convert read1_sequence into one-hot encoding of 4 bits
    # read1_sequence = np.array([[float(nt == x) for x in config.NUC_ORDER] for nt in read1_sequence])
    read1_sequence = np.array([[float(nt != x) for x in config.NUC_ORDER] for nt in read1_sequence])
    
    # add up the counts for each nucleotide
    read1_nt_counts = np.sum(read1_sequence, axis=0)

    # return None if one or more of the nucleotides doesn't show up in the read1_sequence
    if np.min(read1_nt_counts) == 0:
        return None

    # convert intensities to an array and split into read1 and read2 intensities
    read1_intensities = intensities[:config.LEN1 - config.NUM_SKIP - config.NUM_SKIP_2] 
    # read1_intensities = intensities[:config.LEN1 - config.NUM_SKIP]
    read2_intensities = intensities[-1 * config.LEN2:]

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

        return signal#.astype(int)


def get_batch_t_signal(params):
    """
    Calculate normalized t-signals for a batch of sequences.
    Return t-signal values as one long string to be written to a file.
    """

    write_str = ''
    skipped = []

    # iterate through signals in the batch
    for (read_ID, read1, start, signal) in zip(*params):

        signal = signal.reshape((config.LEN1+config.LEN2, 4))

        # impute missing data
        signal = impute(signal, config.NAN_LIMIT)

        # skip if too much missing data
        if signal is None:
            skipped.append(read_ID)
            continue

        # calculate normalization
        normed_signal = get_normalized_intensities(signal, read1)

        # normalize t-signal and return the output as a string
        if normed_signal is not None:
            t_signal = get_t_signal(normed_signal)
            write_str += '{}\t{}\t{}\n'.format(read_ID, str(start),
                                                  '\t'.join(['{:.3f}'.format(x) for x in t_signal]))

        else:
            skipped.append(read_ID)
    
    return skipped, write_str


def calculate_intensities(intensity_file, keep_dict, outdir, num_processes):
    """Iterate through intensity file, extract values for mapping reads, calculate normalized T-signal.
    Writes intensities to file.

    Arguments:
        intensity_file: path to gzipped intensity file
        keep_dict: dictionary of read IDs to sequences
        outdir: directory for output files
        num_processes: number of processes
    """

    # keep track of how many reads we skip due having too many 0's
    skipped = []
    num_reads_kept = 0

    outfile = open(os.path.join(outdir, 'normalized_t_signal.txt'), 'w')

    # if indicated, run parallel version
    if num_processes > 1:
        print("Running parallel version")
        chunks = []
    else:
        print("Running non-parallel version")


    # read intensity file and check if it's in keep_dict
    IDs, read1s, starts, intensity_values = [], [], [], []
    ix = 0
    line_num = 0

    if intensity_file[-3:] == ".gz": infile = gzip.open(intensity_file, 'rt') #is the intensity file .gz?
    elif intensity_file[-4:] == ".txt": infile = open(intensity_file, 'r') #not a gzipped file
    else: raise ValueError("Intensity file does not end in .gz or .txt")
    for line in infile:
        line = line.split()
        line_num += 1
        read_ID = config.intensity_line_to_ID(line)

        try:
            seq, tail_start = keep_dict[read_ID]
            
        except:
            continue

        IDs.append(read_ID)
        read1s.append(seq)
        starts.append(tail_start)
        intensity_values.append(line[config.SIGNAL_COL_START: config.SIGNAL_COL_END])
        ix += 1

        if ix == config.CHUNKSIZE:
            chunk = (IDs, read1s, starts, np.array(intensity_values, dtype=int))
            IDs, read1s, starts, intensity_values = [], [], [], []
            ix = 0

            # if indicated, run parallel version
            if num_processes > 1:
                chunks.append(chunk)

                # once we've read enough chunks, calculate t-signals in parallel
                if len(chunks) == num_processes:
                    full_write_str = ''
                    with concurrent.futures.ProcessPoolExecutor() as executor:
                        results = executor.map(get_batch_t_signal, chunks)

                        # record how many signals were skipped and write results to a file
                        for sk, write_str in results:
                            time.sleep(0.0001)
                            skipped += sk
                            num_reads_kept += write_str.count('\n')
                            full_write_str += write_str

                    outfile.write(full_write_str)

                    full_write_str = ''
                    chunks = []

            # otherwise run sequentially
            else:
                # calculate t-signal
                sk, write_str = get_batch_t_signal(chunk)

                # write to file
                skipped += sk
                num_reads_kept += write_str.count('\n')
                outfile.write(write_str)
    

    # calculate t-signals for the last chunks
    if num_processes > 1:
        if ix > 0:
            chunk = (IDs, read1s, starts, np.array(intensity_values, dtype=int))
            chunks.append(chunk)

        if len(chunks) > 0:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = executor.map(get_batch_t_signal, chunks)
                for sk, write_str in results:
                    skipped += sk
                    num_reads_kept += write_str.count('\n')
                    outfile.write(write_str)
    else:
        if ix > 0:
            chunk = (IDs, read1s, starts, np.array(intensity_values, dtype=int))

            # calculate t-signal
            sk, write_str = get_batch_t_signal(chunk)

            # write to file
            skipped += sk
            num_reads_kept += write_str.count('\n')
            outfile.write(write_str)
    infile.close()
    outfile.close()
    return skipped, num_reads_kept


