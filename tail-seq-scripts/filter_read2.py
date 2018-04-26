from optparse import OptionParser
import regex

import numpy as np
import pandas as pd


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-i", dest="INFILE", help="path to fastq file for read2")
    # parser.add_option("-r", dest="READ1_FILTERED", help="path to filtered information for read1")
    parser.add_option("-m","--manual", dest="MANUAL", help="output file for manually called tails")
    parser.add_option("-o","--outfile", dest="OUTFILE", help="output file")

    (options, args) = parser.parse_args()

    tail_file = open(options.MANUAL, 'w')
    out_file = open(options.OUTFILE, 'w')

    with open(options.INFILE, 'r') as infile:
        while True:
            read_ID = infile.readline()
            if len(read_ID) == 0:
                break

            if read_ID[0] != '@':
                raise ValueError('fastq file entries must start with @')

            read_ID = read_ID[1:].replace('\n','')
            seq = infile.readline()
            _ = infile.readline()
            _ = infile.readline()

            # pass if basecaller confused about more than 2 bases in first 25 nucleotides
            if seq[:25].count('N') >= 2:
                continue

            # look for at least 11 contiguous T's in first 30 nucleotides, allowing 1 error
            match = regex.search(r'TTTTTTTTTTT{e<=1}', seq[:30])

            # if found, keep the ID and where the tail starts
            if match is not None:
                match_start = match.start()
                out_file.write('{}\t{}\n'.format(read_ID, match_start))

            # otherwise manually call tail if read starts with >= 4 contiguous T's
            else:
                if seq[:4] == 'TTTT':
                    TL = 4
                    for nt in seq[4:10]:
                        if nt == 'T':
                            TL += 1
                        else:
                            break
                    tail_file.write('{}\t{}\tmanual\n'.format(read_ID, TL))

    tail_file.close()
    out_file.close()


    # # convert to dataframes
    # keep = pd.DataFrame(keep)
    # keep.columns = ['read_ID', 'tail_start']
    # keep = keep.set_index('read_ID')

    # short_tails = pd.DataFrame(short_tails)
    # short_tails.columns = ['read_ID', 'tail_length']
    # short_tails = short_tails.set_index('read_ID')

    # # import read1 data
    # dtype_dict = {0:str, 1:int, 2:int, 3:str, 4:str, 5:str}
    # read1 = pd.read_csv(options.READ1_FILTERED, sep='\t', header=None, dtype=dtype_dict).set_index(5)
    # keep = pd.concat([read1, keep], axis=1, join='inner')
    # keep.to_csv(options.OUTFILE, sep='\t', header=None)

    # short_tails = pd.concat([read1, short_tails], axis=1, join='inner')
    # short_tails['called'] = 'manual'
    # short_tails.to_csv(options.MANUAL, sep='\t', header=None)

