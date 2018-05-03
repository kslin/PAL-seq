import matplotlib
matplotlib.use('Agg')

from optparse import OptionParser
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-s", dest="INFILE", help="t-signal file")
    parser.add_option("-o", dest="OUTDIR", help="plotting file")
    parser.add_option("-b", dest="BINS", type="int", default=100, help="num bins")

    (options, args) = parser.parse_args()

    # make an empty array for storing the densities
    densities = np.zeros((options.BINS-1, config.LEN2))

    # define the bins
    bin_boundaries = np.linspace(config.LOWERBOUND, config.UPPERBOUND, options.BINS)

    data = pd.read_csv(options.INFILE, sep='\t', header=None)

    # for every t-signal, calculate its density in each bin and add to the array
    for i in range(config.LEN2):
        vals = data[2+i]
        xs, _ = np.histogram(vals, bins=bin_boundaries, normed=True)
        densities[:, i] = xs

    # plot a heatmap
    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    ax.pcolor(densities, cmap=plt.cm.Blues)
    ax.set_yticks(np.linspace(0, options.BINS, 11))
    ax.set_yticklabels(np.linspace(config.LOWERBOUND, config.UPPERBOUND, 11))

    ax.set_ylabel('Normalized T-signal')
    ax.set_xlabel('Position from 3` end')

    fig.savefig(os.path.join(options.OUTDIR, 'signal_plot.pdf'))
    plt.close()
