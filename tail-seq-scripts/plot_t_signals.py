from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config
import helpers


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-f", dest="INFILE", help="t-signal file")
    parser.add_option("-o", dest="OUTFILE", help="plotting file")
    parser.add_option("-b", dest="BINS", type="int", help="num bins") # I used 100
    parser.add_option("-l", dest="LEN", type="int", help="length of t-signal") # e.g. 250

    (options, args) = parser.parse_args()

    densities = np.zeros((options.BINS-1, options.LEN))
    bin_boundaries = np.linspace(config.LOWERBOUND, config.UPPERBOUND, options.BINS)

    data = pd.read_csv(options.INFILE, sep='\t', header=None)

    for i in range(options.LEN):
        vals = data[3+i]
        xs, _ = np.histogram(vals, bins=bin_boundaries, normed=True)
        densities[:, i] = xs

    plt.pcolor(densities, cmap=plt.cm.Blues)
    plt.savefig(options.OUTFILE)
    plt.close()

