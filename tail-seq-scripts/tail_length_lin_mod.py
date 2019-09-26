import concurrent.futures
from optparse import OptionParser
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

import config
import tail_length_helpers

import pdb 

if __name__ == '__main__':

    parser = OptionParser()
    # parser.add_option("-s", dest="SIGNAL", help="signal file output from get_signal.py")
    parser.add_option("-o", "--outdir", dest="OUTDIR", help="output directory")
    # parser.add_option("--twostate", dest="TWOSTATE", default=True, help="toggle for 2-state model")
    # parser.add_option("-f", "--futures", dest="FUTURES", type="int", default=1, help="number of threads to use")

    (options, args) = parser.parse_args()



    print("Writing tail-length outputs to {}".format(options.OUTDIR))

    # find how long the signal file is
    logfile = pd.read_csv(os.path.join(options.OUTDIR, 'logfile.txt'), sep='\t', header=None, index_col=0)
    file_length = int(logfile.loc['Reads for HMM'][1])

    # quit if file too small
    # if file_length < config.TRAIN_SIZE:
    #     print("The training size must be larger than the input file size")
    #     sys.exit()

    # warnings.filterwarnings("ignore", category=DeprecationWarning)
    # warnings.filterwarnings("ignore", category=RuntimeWarning)

    signal_file = os.path.join(options.OUTDIR, 'normalized_t_signal_all.txt')
    signal_file_std = os.path.join(options.OUTDIR, 'normalized_t_signal_stds.txt')
    short_tail_file = os.path.join(options.OUTDIR, 'short_tails.txt')
    all_info_file = os.path.join(options.OUTDIR, 'all_read_info.txt')
    model_params_file = os.path.join(options.OUTDIR, 'model_params.txt')
    std_meds_file = os.path.join(options.OUTDIR, 'standard_medians.txt')

    # if the necesary input files don't exist, quit
    for f in [signal_file, signal_file_std, short_tail_file, all_info_file]:
        if not os.path.exists(f):
            print("{} does not exist. Run T-signal script first with the same output directory.".format(f))
            sys.exit()

    ### get training set ###
    ### Now this training set is just the standards, fitting to a linear model. TJE 2019 09 26###
    t0 = time.time()

    # collect random rows to train
    training_array_dict, std_df = tail_length_helpers.read_training_set(signal_file_std)

    print('{:.3f} seconds'.format(time.time() - t0))
    print("Training Linear Model...")
    t0 = time.time()


    # create linear model, returns an array with the linear model parameters.
    LinRegArr = tail_length_helpers.train_model(training_array_dict,model_params_file,std_meds_file)

    print('{:.3f} seconds'.format(time.time() - t0))
    print("Calculating tail-lengths...")
    t0 = time.time()

    #may not need parallelization.
    all_ids, all_tls = tail_length_helpers.get_batch_tail_length(signal_file,LinRegArr)
    print('{:.3f} seconds'.format(time.time() - t0))


    # creat DataFrame of tail length values
    tail_length_df = pd.DataFrame({'read_ID': all_ids, 'tail_length': all_tls, 'method': 'LM'})
    tail_length_df.loc[tail_length_df['read_ID'].isin(std_df[0]),'method'] = 'LM_STD' #standards called by LM 
    tail_length_df = tail_length_df.set_index('read_ID')

    print('{:d} tail lengths calculated by linear model.'.format(len(tail_length_df)))

    # read in short tails and append
    short_tails = pd.read_csv(short_tail_file, sep='\t', index_col='read_ID')
    short_tails['method'] = 'manual_call'
    tail_length_df = pd.concat([tail_length_df, short_tails])

    # import and merge with other information
    other_info = pd.read_csv(all_info_file, sep='\t', index_col='read_ID')
    merged = pd.concat([tail_length_df, other_info], axis=1, join='inner')

    print('{:d} tail lengths in total.'.format(len(tail_length_df)))

    # write tail lengths to a file
    merged.to_csv(os.path.join(options.OUTDIR, 'tail_lengths.txt'), sep='\t')
