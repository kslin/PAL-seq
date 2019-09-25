import concurrent.futures
from optparse import OptionParser
import os
import sys
import time
import warnings

from hmmlearn import hmm
import numpy as np
import pandas as pd

import config
import tail_length_helpers

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-s", dest="SIGNAL", help="signal file output from get_signal.py")
    parser.add_option("-o", "--outdir", dest="OUTDIR", help="output directory")
    parser.add_option("--twostate", dest="TWOSTATE", default=True, help="toggle for 2-state model")
    parser.add_option("-f", "--futures", dest="FUTURES", type="int", default=1, help="number of threads to use")

    (options, args) = parser.parse_args()
    #parse the boolean here
    if options.TWOSTATE=='True': options.TWOSTATE=True
    elif options.TWOSTATE=='False': options.TWOSTATE=False
    else: options.TWOSTATE=True


    print("Writing tail-length outputs to {}".format(options.OUTDIR))

    # find how long the signal file is
    logfile = pd.read_csv(os.path.join(options.OUTDIR, 'logfile.txt'), sep='\t', header=None, index_col=0)
    file_length = int(logfile.loc['Reads for HMM'][1])

    # quit if file too small
    if file_length < config.TRAIN_SIZE:
        print("The training size must be larger than the input file size")
        sys.exit()

    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    signal_file = os.path.join(options.OUTDIR, 'normalized_t_signal.txt')
    signal_file_std = os.path.join(options.OUTDIR, 'normalized_t_signal.txt')
    short_tail_file = os.path.join(options.OUTDIR, 'short_tails.txt')
    all_info_file = os.path.join(options.OUTDIR, 'all_read_info.txt')

    # if the necesary input files don't exist, quit
    for f in [signal_file, short_tail_file, all_info_file]:
        if not os.path.exists(f):
            print("{} does not exist. Run T-signal script first with the same output directory.".format(f))
            sys.exit()

    ### get training set ###
    t0 = time.time()

    # collect random rows to train
    training_values, training_seq_lengths = tail_length_helpers.read_training_set(signal_file, file_length, config.TRAIN_SIZE, random_seed=0)

    print('{:.3f} seconds'.format(time.time() - t0))
    print("Training HMM...")
    t0 = time.time()

    # create HMM model
    if options.TWOSTATE:
        params = config.HMM_PARAMS_2STATE
        MODEL = hmm.GaussianHMM(n_components=3, init_params="", n_iter=config.MAX_ITER, tol=config.TOL)
        
    else:
        params = config.HMM_PARAMS_3STATE
        MODEL = hmm.GaussianHMM(n_components=4, init_params="", n_iter=config.MAX_ITER, tol=config.TOL)

    MODEL.startprob_ = params['START_PROB_INIT']
    MODEL.transmat_ = params['TRANSITION_INIT']
    MODEL.means_ = params['MEANS_INIT']
    MODEL.covars_ = params['VARS_INIT']

    # train model
    MODEL.fit(training_values, lengths=training_seq_lengths)
    print(MODEL.monitor_)
    print(MODEL.startprob_, MODEL.transmat_, MODEL.means_, MODEL.covars_)

    print('{:.3f} seconds'.format(time.time() - t0))
    print("Calculating tail-lengths...")
    t0 = time.time()

    # define function that calculates tail length
    def get_batch_tail_length(lines):
        """
        Given a chunk of a signal file, use the model to calculate tail lengths
        """
        # extract metadata
        ids, signals, lengths = [], [], []

        # extract signal
        for line in lines:
            line = line[:-1].split('\t') # remove newline character and split by tab
            line = tail_length_helpers.smooth(line)
            ids.append(line[0])
            tail_start = int(float(line[1]))
            signals.append([config.START_SIGNAL])
            for val in line[2 + tail_start:]:
                signals.append([float(val)])
            # if line[0] == '1101:13993:2125':
            #     print("ID found, line 1:")
            #     print(line)
            #     print(line[0])
            #     print(config.START_SIGNAL)
            #     print(tail_start)
            #     print(config.LEN2 - tail_start + 1)
            #     FinalIdentifier = i #Marks the index where the read occurs.
            lengths.append(config.LEN2 - tail_start + 1)

        # predict states
        signals = np.array(signals)
        emissions = MODEL.predict(signals, lengths=lengths)

        # get tail lengths from emissions
        tail_lengths = []
        start_ix = 0
        for length in lengths:
            single_emission = emissions[start_ix: start_ix + length]
            tail_lengths.append(tail_length_helpers.get_tail_length_from_emissions(single_emission, options.TWOSTATE))
            start_ix += length
        #     if i == FinalIdentifier:
        #         print("ID found, line 2:")
        #         print(start_ix)
        #         print(length)
        #         print(single_emission)
        # if '1101:13993:2125' in ids:
        #     print("ID found, line 3:")
        #     print(tail_lengths[ids.index('1101:13993:2125')]) 
        return ids, tail_lengths

    # read in signals as chunks and calculate viterbi emissions
    chunk_iterator = pd.read_csv(signal_file, sep='\t', header=None, iterator=True, chunksize=10000)
    finished = 0
    all_ids, all_tls = [], []

    # if indicated, run parallel version
    if config.FUTURES > 1:
        print("Running parallel version")
    else:
        print("Running non-parallel version")

    ix = 0
    chunks = []
    chunk = []
    with open(signal_file, 'r') as infile:
        for line in infile:

            # add lines to the chunk
            chunk.append(line)
            ix += 1
            finished += 1

            # when we've collected a chunk, calculate the tail lengths
            if ix == config.CHUNKSIZE:

                # if indicated, run in parallel
                if config.FUTURES > 1:
                    chunks.append(chunk)
                    if len(chunks) == config.FUTURES:
                        with concurrent.futures.ProcessPoolExecutor() as executor:
                            results = executor.map(get_batch_tail_length, chunks)

                            # add results to the dataframe
                            for ids, tls in results:
                                all_ids += ids
                                all_tls += tls
                                time.sleep(0.0001)

                        chunks = []

                # otherwise, run sequentially
                else:
                    ids, tls = get_batch_tail_length(chunk)
                    all_ids += ids
                    all_tls += tls
                    break

                ix = 0
                chunk = []

    # process last batch
    if config.FUTURES > 1:
        if len(chunk) > 0:
            chunks.append(chunk)

        if len(chunks) > 0:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = executor.map(get_batch_tail_length, chunks)

                # add results to the dataframe
                for ids, tls in results:
                    all_ids += ids
                    all_tls += tls
                    time.sleep(0.0001)
    else:
        if len(chunk) > 0:
            ids, tls = get_batch_tail_length(chunk)
            all_ids += ids
            all_tls += tls


    # creat DataFrame of tail length values
    tail_length_df = pd.DataFrame({'read_ID': all_ids, 'tail_length': all_tls, 'method': 'HMM'}).set_index('read_ID')

    print(len(tail_length_df))

    # read in short tails and append
    short_tails = pd.read_csv(short_tail_file, sep='\t', index_col='read_ID')
    short_tails['method'] = 'manual_call'
    tail_length_df = pd.concat([tail_length_df, short_tails])

    # import and merge with other information
    other_info = pd.read_csv(all_info_file, sep='\t', index_col='read_ID')
    merged = pd.concat([tail_length_df, other_info], axis=1, join='inner')

    print(len(tail_length_df), len(merged))

    # write tail lengths to a file
    merged.to_csv(os.path.join(options.OUTDIR, 'tail_lengths.txt'), sep='\t')
