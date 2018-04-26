from optparse import OptionParser

import numpy as np
import pandas as pd


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-i", dest="INFILE", help="path to input file with 3p ends")
    parser.add_option("-o","--outfile", dest="OUTFILE", help="output file")
    parser.add_option("-d","--dropped", dest="DROPPED", help="output file for dropped reads")

    (options, args) = parser.parse_args()

    # import data
    data = pd.read_csv(options.INFILE, sep='\t', header=None)

    # split into plus and minus strand reads and sort each
    plus = data[data[3] == '+'].sort_values(1, ascending=False)
    minus = data[data[3] == '-'].sort_values(1, ascending=True)

    # for reads that map to multiple exons of the same transcript, take later exon
    plus = plus.drop_duplicates(subset=[4,5], keep='first')
    minus = minus.drop_duplicates(subset=[4,5], keep='first')

    # discard reads that map to multiple genes
    plus_dedup = plus.drop_duplicates(subset=[5], keep=False)
    minus_dedup = minus.drop_duplicates(subset=[5], keep=False)
    dropped_reads = list(set(list(plus[5])) - set(list(plus_dedup[5]))) + list(set(list(minus[5])) - set(list(minus_dedup[5])))

    print('{} reads dropped due to mapping to multiple genes'.format(len(dropped_reads)))

    # write dropped reads to a file
    with open(options.DROPPED, 'w') as dropped_file:
        dropped_file.write('\n'.join(dropped_reads))

    # concat plus and minus strand information and write to output
    data = pd.concat([plus_dedup, minus_dedup], axis=0)

    data.to_csv(options.OUTFILE, sep='\t', index=False, header=None)
