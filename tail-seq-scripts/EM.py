import concurrent.futures
from optparse import OptionParser
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm, linregress
from sklearn.mixture import GaussianMixture

import helpers

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-s", dest="SIGNAL", help="signal file output from get_signal.py")
    parser.add_option("-l", dest="LEN", type="int", help="length of read2")
    parser.add_option("-t", dest="TRAIN_SIZE", type="int", default=10000, help="size of training set")
    parser.add_option("-m", dest="MAXITER", type="int", default=100, help="maximum number of iterations")
    parser.add_option("--tol", dest="TOL", type="int", default=10, help="tolerance for EM algorithm")
    parser.add_option("-o", "--outfile", dest="OUTFILE", help="output file for tail lengths")
    parser.add_option("-p", "--plot", dest="PLOT", default=None, help="directory to save parameter plots, set as None to skip plotting")
    parser.add_option("-f", "--futures", dest="FUTURES", action="store_true", default=False, help="set to True to run in parallel")

    (options, args) = parser.parse_args()

    if options.FUTURES:
        print("running parallel version")
    else:
        print("running non-parallel version")

    # read in signal data from get_signal.py
    t_signals_df = pd.read_csv(options.SIGNAL,sep='\t').set_index('ID')
    a_mu = np.array(t_signals_df.drop('Tail_start',1).loc['1201:20207:51631'])
    t_signals_df = t_signals_df.iloc[np.random.choice(np.arange(len(t_signals_df)), size=options.TRAIN_SIZE, replace=False)]
    starts = np.array(t_signals_df['Tail_start'])
    t_signals_df = t_signals_df.drop('Tail_start', 1)
    seq_ids = np.array(t_signals_df.index)
    t_signals = t_signals_df.values
    print(t_signals.shape)

    # set initial parameters
    rates = np.linspace(-1,-0.1,250)
    # a_mu = np.array([1.0]*250)
    a_sigma = np.array([0.15]*options.LEN)
    nota_mu1, nota_sigma1 = np.zeros(options.LEN), np.linspace(0.5,0.15,options.LEN)
    nota_mu2, nota_sigma2 = 0, 0.1
    ratio = 0.25

    # empty array for storing tails
    tails = np.zeros(t_signals.shape[0])

    # a_mu, a_sigma, rates = np.array(a_mu_init), a_sigma_init, rates_init
    # nota_mu1, nota_sigma1 = np.zeros(len(a_mu_init)), np.array([0.5]*len(a_mu_init))
    # nota_mu2, nota_sigma2, ratio = 0, 0.1, 0.25

    # plot initial parameters
    if options.PLOT:
        plt.figure(figsize=(15,4))
        plt.scatter(range(len(a_mu)), a_mu, label='a_mu')
        plt.scatter(range(len(a_mu)), a_sigma, label='a_sigma')
        plt.scatter(range(len(nota_mu1)), nota_mu1, label='nota_mu')
        plt.scatter(range(len(a_mu)), nota_sigma1, label='nota_sigma')
        plt.scatter(range(len(rates)), rates, label='rate')
        plt.legend(frameon=False)
        plt.savefig(options.PLOT + 'param_plot_0.pdf')
        plt.close()

    for ix in range(options.MAXITER):
        # create arrays for storing A and non-A signals
        a_vals, nota_vals = [a_mu.reshape(1, options.LEN)], [nota_mu1.reshape(1, options.LEN)]
        rate_vals, new_tails, new_seq_ids = [], [], []

        # if specified, run M-step in parallel
        if options.FUTURES:
            batch_size = 200
            executor = concurrent.futures.ProcessPoolExecutor()
            futures = []
            num_batches = int(np.ceil(len(tails) / batch_size))

            # submit batches of signals to the helper function
            for i in range(num_batches):
                seq_id_batch = seq_ids[i*batch_size:min(len(tails), (i+1)*batch_size)]
                t_signal_batch = t_signals[i*batch_size:min(len(tails), (i+1)*batch_size), :]
                start_batch = starts[i*batch_size:min(len(tails), (i+1)*batch_size)]
                futures.append(executor.submit(helpers.get_tail_length_batch_gmm, seq_id_batch, start_batch, t_signal_batch,
                                               a_mu, a_sigma, nota_mu1, nota_sigma1, nota_mu2, nota_sigma2, ratio, rates))
                time.sleep(0.0001)

            for future in concurrent.futures.as_completed(futures):
                ids, tail_length, A_vals, notA_vals, rate = future.result()
                new_seq_ids += ids
                new_tails += tail_length
                a_vals += A_vals
                nota_vals += notA_vals
                rate_vals += rate

            executor.shutdown()

        # otherwise, run sequentially
        else:
            new_seq_ids = seq_ids
            for i in range(t_signals.shape[0]):
                tail_length, A_vals, notA_vals, rate = helpers.get_tail_length_gmm(t_signals[i,:], starts[i], a_mu, a_sigma,
                                                                               nota_mu1, nota_sigma1, nota_mu2, nota_sigma2, ratio, rates)

                new_tails.append(tail_length)
                a_vals.append(A_vals)
                nota_vals.append(notA_vals)
                rate_vals.append(rate)
        
        # concatenate the results into an array
        new_tails = np.array(new_tails)[np.argsort(new_seq_ids)]
        a_vals = np.concatenate(a_vals, axis=0)
        nota_vals = np.concatenate(nota_vals, axis=0)
        rate_vals = np.concatenate(rate_vals, axis=0)


        # plot a scatter plot of how the tails change from one step to the next
        if options.PLOT is not None:
            plt.figure(figsize=(7,7))
            plt.scatter(tails, new_tails)
            plt.savefig(options.PLOT + 'plot_change_{}.pdf'.format(ix+1))
            plt.close()

        # stop once we reached the tolerance
        diff = np.sum(np.abs(new_tails-tails))
        if diff < options.TOL:
            break

        # otherwise, update the parameters (E-step)
        else:
            tails = new_tails
            a_mu = np.nanmean(a_vals, axis=0)
            nota_mu1 = np.nanmean(nota_vals, axis=0)

            # for i in range(len(nota_mu1)):
            #     a_mu[i] = np.nanmedian(a_vals[:, max(0,i-5): min(len(nota_mu1),i+5)])
            # print(np.nanstd(a_vals, axis=0).shape)
            # a_sigma = np.maximum(np.nanstd(a_vals, axis=0), 0.1)
            a_sigma = np.array([np.nanmean(np.maximum(np.nanstd(a_vals, axis=0), 0.1))]*len(a_mu))
            nota_sigma1 = np.nanstd(nota_vals, axis=0)
            # print(np.nanmedian(np.nanstd(a_vals, axis=0)))

            # adjust first 10 positions of parameters
            a_mu[:20] = [np.nanmean(a_mu[20:30])]*20
            nota_mu1[:20] = [np.nanmean(nota_mu1[20:30])]*20
            nota_sigma1[:20] = [np.nanmean(nota_sigma1[20:30])]*20

            # adjust last 10 positions of parameters
            a_mu[-10:] = [np.nanmean(a_mu[-20:-10])]*10


            ## Uncomment for GMM

            # nota_sigma1_list, nota_mu2_list, nota_sigma2_list, ratio_list = [], [], [], []
            # for i in range(length):
            #     A = GaussianMixture(n_components=2, covariance_type='diag')
            #     vals = [x for x in nota_vals[:, max(0,i-5): min(length, i+5)].flatten() if np.isnan(x) == False]
            #     if len(vals) <= 2:
            #         nota_mu2_list.append(0)
            #         print('no_vals {}'.format(i))
            #         continue
            #     A.fit(np.array(vals).reshape(-1,1))
            #     max_ix = np.argsort(A.means_.reshape(2))
            #     nota_mu1[i] = A.means_[max_ix[1], 0]
            #     nota_sigma1_list.append(A.covariances_[max_ix[1], 0])

            #     nota_mu2_list.append(A.means_[max_ix[0], 0])
            #     nota_sigma2_list.append(A.covariances_[max_ix[0], 0])

            #     ratio_list.append(A.weights_[max_ix[1]])
            
            # nota_sigma1 = np.nanmean(nota_sigma1_list)
            # nota_mu2 = np.nanmean(nota_mu2_list)
            # nota_sigma2 = np.nanmean(nota_sigma2_list)
            # ratio = np.nanmean(ratio_list)
            # print("ratio: {}".format(ratio))
            rates = np.nanmean(rate_vals, axis=0)
            xs, ys = [], []
            for i, r in enumerate(rates):
                if i == 0:
                    continue
                if np.isnan(r):
                    continue
                xs.append(i)
                ys.append(r)
            # rates[0] = np.nanmean(rates)

            if options.PLOT is not None:
                plt.figure(figsize=(15,4))
                plt.scatter(range(len(a_mu)), a_mu, label='a_mu')
                plt.scatter(range(len(a_mu)), a_sigma, label='a_sigma')
                plt.scatter(range(len(a_mu)), nota_mu1, label='nota_mu')
                plt.scatter(range(len(a_mu)), nota_sigma1, label='nota_sigma')
                plt.scatter(xs, ys, label='rates')
                rate_constant, _ = curve_fit(helpers.exp_fit, np.array(xs), np.array(ys))
                rates = np.array([-1] + list(np.exp(-1*rate_constant[0] / np.arange(1,options.LEN)) - 1))
                plt.scatter(range(len(rates)), rates, label='rates_curve')
                plt.title('{}'.format(diff))
                plt.legend(frameon=False)
                plt.savefig(options.PLOT + 'param_plot_{}.pdf'.format(ix+1))
                plt.close()

    # write tail lengths to a file
    tail_df = pd.DataFrame({'ID': sorted(new_seq_ids), 'Tail_length': new_tails})
    tail_df.to_csv(options.OUTFILE, sep='\t', index=False)

    with open(options.PLOT + 'params.txt','w') as f:
        f.write('\t'.join([str(x) for x in a_mu]) + '\n')
        f.write('\t'.join([str(x) for x in a_sigma]) + '\n')
        f.write('\t'.join([str(x) for x in nota_mu1]) + '\n')
        f.write('\t'.join([str(x) for x in nota_sigma1]) + '\n')
        f.write(str(rate_constant))
        
