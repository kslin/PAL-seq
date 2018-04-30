import gzip

import numpy as np
import pandas as pd
import regex

import config


def reverse_complement(seq):
    """Get reverse complement of sequence"""
    nt_dict = {'A':'T', 'T': 'A', 'C': 'G', 'G':'C', 'N': 'N'}
    return ''.join([nt_dict[nt] for nt in seq][::-1])


def dedup_bed(bedfile):
    """Deduplicate intersectBed results.
    If a read maps to multiple genes, discard.
    If a read maps to multiple exons of a gene, use later exon.
    Extract 3' end of read. Note that reads are reverse complement of genes.

    Arguments:
        bedfile: pandas DataFrame of intersectBed columns [chr, start, end, read_ID, strand, gene]
    """

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
    """Parse fastq1 file, extract reads that match standard sequences.

    Arguments:
        fastq1: path to gzipped fastq file for read1
        keep_dict: dictionary of {read_ID: None} for reads that pass filters so far
        standard_dict: dictionary of standard sequences
    """
    new_keep_dict = {}
    standard_reads = []
    with gzip.open(fastq1, 'rt') as infile:
        while True:
            line1 = infile.readline()
            if line1 == '':
                break

            if line1[0] != '@':
                raise ValueError('fastq file entries must start with @')

            read_ID = config.fastq_header_to_ID(line1[1:]) # exclude leading @
            seq = infile.readline()[:-1] # exclude newline character
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
    """Parse fastq2 file, filter low quality or very short tails.
    Manually call tails between 4 and 9 nucleotides (inclusive).

    Arguments:
        fastq2: path to gzipped fastq file for read1
        keep_dict: dictionary of {read_ID: read1 sequence} for reads that pass filters so far
        short_tail_outfile_path: path to output file for manually-called short tails
    """
    new_keep_dict = {}
    short_tail_outfile = open(short_tail_outfile_path, 'w')
    dropped_read2 = []
    with gzip.open(fastq2, 'rt') as infile:
        while True:
            line1 = infile.readline()
            if line1 == '':
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

