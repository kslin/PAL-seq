import matplotlib
matplotlib.use('Agg')

from optparse import OptionParser
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config


def get_cdf(vals):
    if len(vals) < 5:
        return [], []

    num_bins = int(len(vals)/5)
    counts, bin_edges = np.histogram(vals, bins=num_bins)
    counts = counts / float(sum(counts))
    return bin_edges[1:], np.cumsum(counts)


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-i", dest="INFILE", help="tail length file")
    parser.add_option("-o", dest="OUTDIR", help="plotting file")

    (options, args) = parser.parse_args()

    # tails = pd.read_csv(options.INFILE, sep='\t', header=None)

    # print(tails.head())

    # standard_names = ['NM_A10_6','NM_A30_6','NM_A110_6','NM_A210_6','NM_A10_7','NM_A50_7','NM_A100_7','NM_A150_7','NM_A200_7','NM_A250_7','NM_A300_7']

    # fig = plt.figure()
    # ax = plt.subplot(1,1,1)

    # standards = tails[tails[0].isin(standard_names)]
    # print(standards.head())
    # for standard, group in standards.groupby(0):
    #     bins, cdf = get_cdf(group[3].values)
    #     ax.plot(bins, cdf, label=standard.replace('NM_',''))

    # ax.legend(loc='lower right', fontsize=7, ncol=2)

    # fig.savefig(os.path.join(options.OUTDIR, 'standard_plot.pdf'))
    # plt.close()

    tails = pd.read_csv(options.INFILE, sep='\t')

    gene_summary = []
    for accession, group in tails.groupby('accession'):
        gene_summary.append([accession, len(group), np.mean(group['tail_length']), np.median(group['tail_length'])])

    gene_summary = pd.DataFrame(gene_summary)
    gene_summary.columns = ['accession', 'num_tags', 'mean_tail_length', 'median_tail_length']
    gene_summary.to_csv(os.path.join(options.OUTDIR, 'tail_lengths_summarized.txt'), sep='\t', index=False)

    fig = plt.figure()
    ax = plt.subplot(1,1,1)

    standards = tails[tails['chr'] == 'standard']
    for standard, group in standards.groupby('accession'):
        bins, cdf = get_cdf(group['tail_length'].values)
        ax.plot(bins, cdf, label=standard)

    ax.legend(loc='lower right', fontsize=7, ncol=2)

    fig.savefig(os.path.join(options.OUTDIR, 'standard_plot.pdf'))
    plt.close()



