import gzip
import regex

import numpy as np
import pandas as pd

import config


def dedup_bed(bedfile):

    # split into plus and minus strand reads and sort each by 3' end
    columns = ['chr','3p_end','read_ID','strand','gene']
    plus = bedfile[bedfile[18] == '+'][[0,2,3,18,20]].sort_values(2, ascending=False)
    plus.columns = columns
    minus = bedfile[bedfile[18] == '-'][[0,1,3,18,20]].sort_values(1, ascending=True)
    minus.columns = columns

    # for reads that map to multiple exons of the same transcript, take later exon
    plus_last_exon = plus.drop_duplicates(subset=['read_ID','gene'], keep='first')
    minus_last_exon = minus.drop_duplicates(subset=['read_ID','gene'], keep='first')

    # discard reads that map to multiple genes
    plus = plus.drop_duplicates(subset=['read_ID'], keep=False)
    minus = minus.drop_duplicates(subset=['read_ID'], keep=False)

    # concat plus and minus strand information
    reads_dedup = pd.concat([plus, minus], axis=0)

    # get dropped reads
    dropped_reads = list(set(bedfile[3].values) - set(reads_dedup['read_ID'].values))

    return reads_dedup, dropped_reads


def parse_read1(fastq1, keep_dict, standard_dict):
    new_keep_dict = {}
    standard_reads = []
    with gzip.open(fastq1, 'r') as infile:
        while True:
            line1 = infile.readline()
            if len(line1) == 0:
                break

            if line1[0] != '@':
                raise ValueError('fastq file entries must start with @')

            read_ID = config.fastq_header_to_ID(line1[1:])
            seq = infile.readline().replace('\n','')
            _ = infile.readline()
            _ = infile.readline()

            is_standard = False
            for standard, name in standard_dict.items():
                if standard in seq:
                    standard_reads.append(['standard', 0, read_ID, '+', standard_dict[standard]])
                    new_keep_dict[read_ID] = seq
                    is_standard = True
                    break

            if is_standard == False:
                if read_ID in keep_dict:
                    new_keep_dict[read_ID] = seq

    if len(standard_reads) == 0:
        standard_reads = pd.DataFrame(None)
    else:
        standard_reads = pd.DataFrame(standard_reads)
        standard_reads.columns = ['chr','3p_end','read_ID','strand','gene']
        standard_reads = standard_reads

    return new_keep_dict, standard_reads


def parse_read2(fastq2, keep_dict, short_tail_outfile_path):
    new_keep_dict = {}
    short_tail_outfile = open(short_tail_outfile_path, 'w')
    dropped_read2 = []
    with gzip.open(fastq2, 'r') as infile:
        while True:
            line1 = infile.readline()
            if len(line1) == 0:
                break

            if line1[0] != '@':
                raise ValueError('fastq file entries must start with @')

            read_ID = config.fastq_header_to_ID(line1[1:])
            seq = infile.readline()
            _ = infile.readline()
            _ = infile.readline()

            if read_ID in keep_dict:

                # pass if basecaller confused about more than 2 bases in first 25 nucleotides
                if seq[:25].count('N') >= 2:
                    dropped_read2.append([read_ID, 'low_qual_read2'])
                    continue

                # look for at least 11 contiguous T's in first 30 nucleotides, allowing 1 error
                match = regex.search(r'TTTTTTTTTTT{e<=1}', seq[:30])

                # if found, keep the ID and where the tail starts
                if match is not None:
                    match_start = match.start()
                    new_keep_dict[read_ID] = (keep_dict[read_ID], match_start)

                # otherwise manually call tail if read starts with >= 4 contiguous T's
                else:
                    if seq[:4] == 'TTTT':
                        TL = 4
                        for nt in seq[4:10]:
                            if nt == 'T':
                                TL += 1
                            else:
                                break
                        short_tail_outfile.write('{}\t{}\n'.format(read_ID, TL))

                    else:
                        dropped_read2.append([read_ID, 'no_tail'])

    return new_keep_dict, dropped_read2

