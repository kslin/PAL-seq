import matplotlib
matplotlib.use('Agg')

from optparse import OptionParser
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config


def get_cdf(vals):
    ### Calculate CDF values for plotting###
    if len(vals) < 5:
        return [], []

    num_bins = int(len(vals)/5)
    counts, bin_edges = np.histogram(vals, bins=num_bins)
    counts = counts / float(sum(counts))

    return bin_edges[1:], np.cumsum(counts)


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-o", dest="OUTDIR", help="plotting file")
    parser.add_option("-f", dest="FILTER", type = 'int', default = 0)

    (options, args) = parser.parse_args()

    # if the necesary input file doesn't exist, quit
    infile = os.path.join(options.OUTDIR, 'tail_lengths.txt')
    if not os.path.exists(infile):
        print("{} does not exist. Run tail_length_hmm.py first with the same output directory.".format(infile))
        sys.exit()

    # read in tail lengths and aggregate by accession
    tails = pd.read_csv(infile, sep='\t')

    gene_summary = []

    #filtering tails by min length in options. 
    tails = tails[tails['tail_length'] >= options.FILTER]

    for accession, group in tails.groupby('accession'):
        gene_summary.append([accession, len(group), np.mean(group['tail_length']), np.median(group['tail_length'])])

    gene_summary = pd.DataFrame(gene_summary)
    gene_summary.columns = ['accession', 'num_tags', 'mean_tail_length', 'median_tail_length']
    gene_summary.to_csv(os.path.join(options.OUTDIR, 'tail_lengths_summarized.txt'), sep='\t', index=False)

    # plot standards if they exist
    standards = tails[tails['chr'] == 'standard']
    if len(standards) > 0:
        fig = plt.figure()
        ax = plt.subplot(1,1,1)

        for standard, group in standards.groupby('accession'):
            bins, cdf = get_cdf(group['tail_length'].values)
            ax.plot(bins, cdf, label=standard)

        ax.legend(loc='lower right', fontsize=7, ncol=2)

        fig.savefig(os.path.join(options.OUTDIR, 'standard_plot.png'))
        plt.close()
