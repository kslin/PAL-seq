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
    parser.add_option("-l", dest="LEN", type="int", help="length of read2")
    parser.add_option("-o", "--outdir", dest="OUTDIR", help="output directory")
    parser.add_option("-t", dest="TRAIN_SIZE", type="int", default=10000, help="size of training set")
    parser.add_option("-m", dest="MAXITER", type="int", default=10000, help="maximum number of iterations")
    parser.add_option("--tol", dest="TOL", type="float", default=0.01, help="tolerance for EM algorithm")
    parser.add_option("-f", "--futures", dest="FUTURES", type="int", default=0, help="number of threads to use")

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
    t_signals_train.columns = ['ID','Gene_name','Tail_start'] + [str(i) for i in range(options.LEN)]
    t_signals_train = t_signals_train.set_index('ID')

    # read in signal data from get_signal.py
    starts = np.array(t_signals_train['Tail_start'])
    t_signals_train = t_signals_train.drop(['Gene_name','Tail_start'], 1)
    seq_ids = np.array(t_signals_train.index)
    t_signals_train = t_signals_train.values.tolist()

    # truncate signal to where the T's start and add the starting state value
    t_signals_train = [[config.START_SIGNAL] + x[s:] for (x,s) in zip(t_signals_train, starts)]

    print("Training HMM")

    

    # create HMM model
    F = ghmm.Float()

    MODEL = ghmm.HMMFromMatrices(F, ghmm.GaussianMixtureDistribution(F),
                                 config.TRANSITIONMATRIX, config.EMISSIONMATRIX, config.PI)

    # MODEL = ghmm.HMMFromMatrices(F, ghmm.GaussianDistribution(F),
    #                              config.TRANSITIONMATRIX2, config.EMISSIONMATRIX2, config.PI2)

    trainingset = ghmm.SequenceSet(F, t_signals_train)
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
        signals = chunk[range(3, options.LEN + 3)].values.tolist()
        signals = [[config.START_SIGNAL] + x[s:] for (x,s) in zip(signals, starts)]
        signals = ghmm.SequenceSet(F, signals)
        emissions = MODEL.viterbi(signals)[0]

        # get tail lengths from HMM emissions
        tail_lengths = [helpers.get_tail_length_from_emissions(x) for x in emissions]
        tail_length_df = pd.DataFrame({'ID': ids, 'Gene_name': genes, 'Start': starts, 'Tail_length': tail_lengths})

        return tail_length_df

    # print(MODEL)
    # help(MODEL)

    # test_signal_list = [2.49865462652, 2.5388379017, 2.27933695711, 1.72512733435, 1.61263282697, 1.93170801239, 1.28087088083, 1.48817647657, 1.34042464698, 1.08774160464, -1.14193026487, -5.0, -2.41408852714, 0.57970458899, 1.65722444695, 1.47895470551, -0.578865950806, -4.25286196628, -2.16752100269, -5.0, -4.69230908126, -5.0, -5.0, -2.01361603913, -1.65742634481, -5.0, -5.0, -5.0, 0.190077227944, 1.15255763343, 0.998999801903, -3.69076254452, -2.16977500495, -4.55282508627, -5.0, -2.73450813159, -0.182450123945, -2.78114329713, -5.0, -0.502110901489, -0.631052066972, -2.8818832548, -2.11142400424, -2.42975589587, -0.933568978346, -2.78100418613, -2.45619865077, -0.34442004885, -2.49496835192, -1.60634276061, -3.60617897949, -3.17591254728, -0.595721375705, 0.243369030871, -2.24822341188, -0.427021100616, 0.40297749154, -0.725096443945, 0.312199184966, 0.562986522956, -1.50045244007, -0.0633483077502, 0.618673124089, 0.798398108669, 0.905505164008, 0.937183074734, 1.25054084553, -0.201025687023, -1.66886040061, -2.25254945516, -2.94075186574, -3.5194776366, -3.05491529893, -1.1036448099, -1.46457232423, -1.91510613575, -2.18525137932, -1.24775317619, -1.99080569581, -2.60639844032, -3.85177778381, -3.96309349952, -1.66021288366, -2.44269815624, -3.29194207266, -2.66436189635, -2.49792777578, -3.07784208667, -1.30828038281, -0.253183661774, 0.320823843292, -0.652881403682, -2.05234265648, -1.09246783732, -1.82022382994, -1.11687057876, -1.35054440246, -0.742606368481, 0.0440101233668, 0.555526857149, 0.174474014989, -0.986877005432, -2.12960167612, -3.35041059577, -2.32199661851, -1.13765367149, -0.273760755103, 0.46772022426, -0.749436278621, -1.52495769884, -1.35866164361, -1.74715360371, -1.89648173594, -1.84462390466, -1.11534306533, -0.571669517198, -0.426285856513, -1.21872848999, -1.81422835645, -1.05324004276, -1.2305591091, -0.994056376697, -0.328523632914, -1.35327739081, -0.78028125844, -0.539252046977, -0.919084377529, -1.34474244078, -1.60014389533, -2.20547685054, -2.49606522235, -3.27950881101, -4.23478471593, -3.14842909308, -1.7585528976, -1.05332166326, -1.65974296245, -2.67082260488, -2.35757174132, -2.141664464, -1.85093725393, -1.60879786602, -1.36548185371, -1.50705819996, -1.83894761627, -1.95486509746, -1.66614048728, -1.9796184686, -2.41810694136, -2.38613684142, -1.70228202893, -2.17225726147, -2.02450829118, -1.45634874245, -1.60681904165, -1.61090264138, -1.60167246679, -2.13798896692, -1.90153497442, -2.39259154785, -2.95030913071, -2.19924387771, -1.81985038299, -1.71619039428, -1.74077857173, -1.45136416821, -1.50839209045, -1.58155751887, -0.560203127551, -0.801782363183, -1.60218426256, -1.81639158861, -1.734660825, -1.31701203563, -2.37175932451, -1.85184822991, -2.48774938961, -2.05529551369, -2.2789085329, -2.2356557482, -1.69906976959, -0.736936596098, -1.35988009852, -2.07942785927, -1.89385215828, -1.59771371596, -1.88357461566, -1.05710008584, -1.40354602227, -1.35988776556, -1.70648261277, -1.12918382405, -1.19277223583, -1.3753973206, -1.21546438356, -1.42612941536, -1.15972736311, -1.42802960509, -1.85830636652, -1.46471547999, -1.25428531688, -0.789300364863, -1.08980173559, -1.12173590569, -1.10862965362, -1.16317983237, -1.34510916348, -1.02540570629, -1.29610690411, -1.26910402721, -1.32451380484, -1.03517718799, -1.99263150292, -2.51635024751, -2.21705364136, -3.01699131846, -2.1790236219, -3.29019899406, -2.58662516066, -2.71146013242, -2.52709017307, -1.76001255143, -1.72230429752, -2.45867749544, -1.6848593148, -3.01910405485, -2.81806600919, -2.4818379897, -2.95584550132, -2.57804812829, -1.97058099538, -1.13381413391, -2.67290149692, -2.40340403313, -2.61427253661, -2.65160946244, -2.87471906002, -1.5972966831, -0.76422534543, -2.29891017459, -1.22914244966, -1.94537535101, -1.99326354198, -2.78037746692, -1.43022980079, -1.30975845746, -1.36524227245, -1.51390536472, -1.02916183082, -2.77300323333]
    # test_signal_list = [[config.START_SIGNAL] + test_signal_list]
    # expected = [4] + [0]*10 + [3]*240
    
    # test_signal = ghmm.SequenceSet(F, test_signal_list)
    # emissions = MODEL.viterbi(test_signal)[0]
    # print(emissions)
    
    # # logprob1, logprob2 = 0, 0
    # # for i in range(30):
    # #     logprob1 = np.log(MODEL.getEmissionProbability(test_signal_list[0][i], emissions[i]))
    # #     logprob2 = np.log(MODEL.getEmissionProbability(test_signal_list[0][i], expected[i]))

    # #     print(test_signal_list[0][i], emissions[i], expected[i], logprob1, logprob2)
    # # print(MODEL.pathPosterior(expected, test_signal))

    # sys.exit()


    print("Calculating tail-lengths")

    # create dataframe to store tail length information
    tail_length_df = pd.DataFrame({'ID':[], 'Gene_name': [], 'Start': [], 'Tail_length': []})

    # read in signals as chunks and calculate viterbi emissions
    chunk_iterator = pd.read_csv(options.SIGNAL, sep='\t', header=None, iterator=True, chunksize=10000)
    finished = 0

    # if indicated, run parallel version
    if options.FUTURES > 0:
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


    else:

        print("Running non-parallel version")
        
        for chunk in chunk_iterator:
            temp = get_batch_tail_length(chunk)
            tail_length_df = pd.concat([tail_length_df, temp])

            finished += len(chunk)
            print('Calculated {} out of {}'.format(finished, file_length))
            # break

    tail_length_df.to_csv(os.path.join(options.OUTDIR, 'tail_lengths_ghmm2.txt'), sep='\t', index=False)

    # get median tail-lengths and write to a separate file
    tail_length_df['Num tags'] = [1]*len(tail_length_df)
    tail_length_df = tail_length_df.groupby('Gene_name').agg({'Tail_length': [np.nanmedian, np.nanmean], 'Num tags': len})
    tail_length_df.to_csv(os.path.join(options.OUTDIR, 'median_tail_lengths_ghmm2.txt'), sep='\t')
