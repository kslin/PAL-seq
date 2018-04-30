import concurrent.futures
from optparse import OptionParser
import os
import sys
import time

import ghmm
import numpy as np
import pandas as pd

import config
import helpers


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-s", dest="SIGNAL", help="signal file output from get_signal.py")
    parser.add_option("-o", "--outdir", dest="OUTDIR", help="output directory")
    parser.add_option("--twostate", action="store_true", dest="TWOSTATE", default=False, help="toggle for 2-state model")
    parser.add_option("-t", dest="TRAIN_SIZE", type="int", default=10000, help="size of training set")
    parser.add_option("-m", dest="MAXITER", type="int", default=10000, help="maximum number of iterations")
    parser.add_option("--tol", dest="TOL", type="float", default=0.01, help="tolerance for EM algorithm")
    parser.add_option("-f", "--futures", dest="FUTURES", type="int", default=1, help="number of threads to use")

    (options, args) = parser.parse_args()

    if options.OUTDIR[-1] == '/':
        options.OUTDIR = options.OUTDIR[:-1]

    print("Writing tail-length outputs to {}".format(options.OUTDIR))

    # find how long the signal file is
    file_length = 0
    with open(options.SIGNAL,'r') as f:
        for line in f:
            file_length += 1

    # quit if file too small
    if file_length < options.TRAIN_SIZE:
        print("The training size must be larger than the input file size")
        sys.exit()

    # if the output directory doesn't exist, quit
    if not os.path.exists(options.OUTDIR):
        print("{} does not exist. Run make-signal with the same output directory".format(options.OUTDIR))
        sys.exit()

    ### get training set ###
    print("Reading in training set")

    # generate which rows to skip when reading in data
    np.random.seed(0)
    skip_ix = np.sort(np.random.choice(range(file_length), size=(file_length - options.TRAIN_SIZE), replace=False))

    # read in data, skipping all rows except the training rows
    t_signals_train = pd.read_csv(options.SIGNAL, sep='\t', skiprows=skip_ix)#.set_index('ID')
    t_signals_train.columns = ['ID','Tail_start'] + [str(i) for i in range(config.LEN2)]
    t_signals_train = t_signals_train.set_index('ID')

    # read in signal data from get_signal.py
    starts = t_signals_train['Tail_start'].values
    t_signals_train = t_signals_train.drop(['Tail_start'], 1)
    seq_ids = np.array(t_signals_train.index)
    t_signals_train = t_signals_train.values.tolist()

    # truncate signal to where the T's start and add the starting state value
    t_signals_train = [[config.START_SIGNAL] + x[s:] for (x,s) in zip(t_signals_train, starts)]

    print("Training HMM")

    # create HMM model
    F = ghmm.Float()

    if options.TWOSTATE:
        params = config.HMM_PARAMS_2STATE
    else:
        params = config.HMM_PARAMS_3STATE

    MODEL = ghmm.HMMFromMatrices(F, ghmm.GaussianMixtureDistribution(F),
                                     params["TRANSITIONMATRIX"], params["EMISSIONMATRIX"], params["PI"])


    trainingset = ghmm.SequenceSet(F, t_signals_train)

    # train model
    MODEL.baumWelch(trainingset, options.MAXITER, options.TOL)

    def get_batch_tail_length(chunk):
        """
        Given a chunk of a signal file, use the model to calculate tail lengths
        and return the results as a data frame
        """
        # extract metadata
        ids = list(chunk[0])
        genes = list(chunk[1])
        starts = list(chunk[2])

        # extract signal
        signals = chunk[np.arange(3, options.LEN + 3)].values.tolist()
        signals = [[config.START_SIGNAL] + x[s:] for (x,s) in zip(signals, starts)]
        signals = ghmm.SequenceSet(F, signals)
        emissions = MODEL.viterbi(signals)[0]

        # get tail lengths from HMM emissions
        tail_lengths = [helpers.get_tail_length_from_emissions(x) for x in emissions]
        tail_length_df = pd.DataFrame({'ID': ids, 'Gene_name': genes, 'Start': starts, 'Tail_length': tail_lengths})

        return tail_length_df

    print("Calculating tail-lengths")

    # create dataframe to store tail length information
    tail_length_df = pd.DataFrame({'ID':[], 'Gene_name': [], 'Start': [], 'Tail_length': []})

    # read in signals as chunks and calculate viterbi emissions
    chunk_iterator = pd.read_csv(options.SIGNAL, sep='\t', header=None, iterator=True, chunksize=10000)
    finished = 0

    # if indicated, run parallel version
    if options.FUTURES > 1:
        print("Running parallel version")

        chunks = []
        t0 = time.time()

        for chunk in chunk_iterator:
            chunks.append(chunk)
            finished += len(chunk)

            if len(chunks) == options.FUTURES:
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    results = executor.map(get_batch_tail_length, chunks)

                    # add results to the dataframe
                    for temp in results:
                        tail_length_df = pd.concat([tail_length_df, temp])
                        time.sleep(0.0001)

                print('Calculated {} out of {}, ({:.2} sec)'.format(finished, file_length, time.time()-t0))
                t0 = time.time()

                chunks = []

        # process last batch
        if len(chunks) > 0:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = executor.map(get_batch_tail_length, chunks)

                # add results to the dataframe
                for temp in results:
                    tail_length_df = pd.concat([tail_length_df, temp])
                    time.sleep(0.0001)

    # otherwise, run sequentially
    else:
        print("Running non-parallel version")
        
        for chunk in chunk_iterator:
            temp = get_batch_tail_length(chunk)
            tail_length_df = pd.concat([tail_length_df, temp])

            finished += len(chunk)
            print('Calculated {} out of {}'.format(finished, file_length))

    # write tail lengths to a file
    tail_length_df.to_csv(os.path.join(options.OUTDIR, 'tail_lengths.txt'), sep='\t', index=False)

    # get median tail lengths and write to a separate file
    tail_length_df['Num tags'] = [1]*len(tail_length_df)
    tail_length_df = tail_length_df.groupby('Gene_name').agg({'Tail_length': [np.median, np.mean], 'Num tags': len})
    tail_length_df.to_csv(os.path.join(options.OUTDIR, 'median_tail_lengths.txt'), sep='\t')
