import gzip
import os

import numpy as np
import pandas as pd
import regex

import config
from itertools import groupby
import pysam
import pdb

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
    columns = ['chr','3p_end','read_ID','strand','accession']
    plus = bedfile[bedfile[18] == '+'][[0,2,3,18,20]].sort_values(2, ascending=False)
    plus.columns = columns
    minus = bedfile[bedfile[18] == '-'][[0,1,3,18,20]].sort_values(1, ascending=True)
    minus.columns = columns

    # for reads that map to multiple exons of the same transcript, take later exon
    plus = plus.drop_duplicates(subset=['read_ID','accession'], keep='first') #Previously, I think these reads were lost. TJE 2019 03 13
    minus = minus.drop_duplicates(subset=['read_ID','accession'], keep='first')

    # discard reads that map to multiple genes
    plus = plus.drop_duplicates(subset=['read_ID'], keep=False)
    minus = minus.drop_duplicates(subset=['read_ID'], keep=False)

    # concat plus and minus strand information
    reads_dedup = pd.concat([plus, minus], axis=0)

    # get dropped reads
    dropped_reads = list(set(bedfile[3].values) - set(reads_dedup['read_ID'].values))

    return reads_dedup, dropped_reads

def parse_read2_BAM(outdir):
    """Parse a BAM file for read2, returning a dictionary of soft clipping.

    Arguments:
        unopened outdir
    """
    sam = pysam.AlignmentFile(os.path.join(outdir, 'Read2STAR_Aligned.sortedByCoord.out.bam'), "rb")
    # sam = pysam.AlignmentFile(os.path.join(outdir, 'A10ShortTailA100.bam'), "rb")
    softClippingDict = {}
    for read in sam:
        read_ID = config.fastq_header_to_ID(read.query_name)
        if read.flag & 16 == 16: 
            # seqRaw = reverse_complement(read.seq) #extract the strand from the bitwise operator.
            softClipping = read.cigartuples[-1][1] #cigar identifiers for soft clipping. 
        else: 
            # seqRaw =  read.seq
            softClipping = read.cigartuples[0][1]
        softClippingDict[read_ID] = softClipping

    return(softClippingDict)

def parse_read1(fastq1, keep_dict, standard_dict):
    """Parse fastq1 file, extract reads that match standard sequences.

    Arguments:
        fastq1: openened fastq1 file, from tarfile, with lines as binary. 
        keep_dict: dictionary of {read_ID: None} for reads that pass filters so far
        standard_dict: dictionary of standard sequences
    """
    new_keep_dict = {}
    standard_reads = []
    line_counter = 0
    for line in fastq1:
        line = line.decode("utf-8") #This line is to decode the tarfile output to utf-8. 
        if line_counter == 0:
            read_ID = config.fastq_header_to_ID(line[1:])
        elif line_counter == 1:
            seq = line[:-1]

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

        line_counter = (line_counter + 1) % 4

    if len(standard_reads) == 0:
        standard_reads = pd.DataFrame(None)
    else:
        standard_reads = pd.DataFrame(standard_reads)
        standard_reads.columns = ['chr','3p_end','read_ID','strand','accession']
        standard_reads = standard_reads

    return new_keep_dict, standard_reads


def parse_read2(fastq2, keep_dict, outdir, softClippingDict, standard_reads, qual_filter = True):
    """Parse fastq2 file, filter low quality or very short tails.
    Manually call tails between 4 and 9 nucleotides (inclusive).

    Arguments:
        fastq2: openened fastq2 file, from tarfile, with lines as binary. 
        keep_dict: dictionary of {read_ID: read1 sequence} for reads that pass filters so far
        outdir: path to output directory

    Modified by TJE on 20190125 to allow argument to turn off low quality filtering. 
    """

    new_keep_dict = {}
    short_tail_outfile = open(os.path.join(outdir, 'short_tails.txt'), 'w')
    short_tail_outfile.write('read_ID\ttail_length\n')
    num_short_tails = 0
    dropped_read2 = []
    line_counter = 0
    standards = set(standard_reads['read_ID'].unique())
    total_lines = 0
    strMove = 30 #What is the thershold for going into the HMM?

    # if qual_filter and config.TRIM_BASES != 4: 
        # r = regex.compile('(%s){e<=1}' % 'TTTTTTTTTTT')

    # elif not qual_filter and config.TRIM_BASES != 4: 
        # r = regex.compile('(%s){e<=2}' % 'TTTTTTTTTTT')

    # elif config.TRIM_BASES == 4:
        # SHOULD BE: If there are >= 30 As and Ts, 
        # r = regex.compile('(%s){e<=1}' % ('T'*strMove)) #30 for this. 

    for line in fastq2:
        line = line.decode("utf-8")
        total_lines += 1
        # if total_lines % 1000 == 0: print("lines processed: %s"%(total_lines))
        if line_counter == 0:
            if line[0] != '@':
                raise ValueError('Fastq2 file must begin with @')
            read_ID = config.fastq_header_to_ID(line[1:])
        elif line_counter == 1:
            seq = line[config.TRIM_BASES:] #In a splint run, this is 0. In a direct lig run, this is 4. 
        line_counter = (line_counter + 1) % 4
        if line_counter != 0 or read_ID not in keep_dict: continue
        # pass if basecaller confused about more than 2 bases in first 25 nucleotides
        if qual_filter and seq[:25].count('N') >= 2:
            dropped_read2.append([read_ID, 'low_qual_read2'])

        else:
            for i in range(0,11):
                r = regex.compile('(%s){e<=1}' % ('A'*(i) + 'T'*(strMove - i))) #3
                # look for at least 11 contiguous T's in first 30 nucleotides, allowing 1 er
                match = r.search(seq[:30]) #changes this to two mismatches if low_qual allow
                # if found, keep the ID and where the tail starts
                if match is not None:
                    match_start = match.start() + i #CHANGED TJE 2019 06 14
                    if seq[match_start] != 'T': match_start += 1 #if error is first pos.
                    if seq[match_start] != 'T': print("Error with starting position.")
                    match_start = match_start + config.TRIM_BASES #CHANGED TJE 2019 06 14
                    new_keep_dict[read_ID] = (keep_dict[read_ID], match_start)
                    break #none of the next statements will be evaluated.
            else: 
                if read_ID in softClippingDict:
                    #All of these reads are in the Read2STAR aligned file.
                    # if read_ID == '1207:14335:91178': pdb.set_trace()

                    #### These next two lines used to remove things that ended in tails if they matched the genome.
                    # if softClippingDict[read_ID] <= 4: #If no soft clipping, no tail.
                    #     dropped_read2.append([read_ID, 'no_tail'])
                    # else: #There is soft clipping and read2 is mapped. 
                    #     #List of tuples nt then count
                    firstBasesList = [[k, len(list(g))] for k, g in groupby(seq)]
                    if firstBasesList[0][0] == 'T':
                        TL = firstBasesList[0][1]
                        short_tail_outfile.write('{}\t{}\n'.format(read_ID, TL))
                        num_short_tails += 1
                    elif len(firstBasesList) > 1 and firstBasesList[0][0] == 'A' and firstBasesList[1][0] == 'T': #uridylation
                        TL = firstBasesList[1][1] #Just the number of Ts
                        short_tail_outfile.write('{}\t{}\n'.format(read_ID, TL))
                        num_short_tails += 1
                    elif len(firstBasesList) > 1 and firstBasesList[0][1] <= 2 and firstBasesList[1][0] == 'T': #allow guanylation for fewer than 2 Gs.
                        TL = firstBasesList[1][1] #Just the number of Ts
                        short_tail_outfile.write('{}\t{}\n'.format(read_ID, TL))
                        num_short_tails += 1
                    #look for tails >= 8 and <19 nt anywhere, no mismatches, must be soft clipped. 
                    else:
                        for i in range(19,7,-1):
                            query = ('T'*i)
                            #This next line requires that the read needs to be soft clipped by at least 3 nt (a 5 nt untemplated tail).
                            if query in seq[:30] and softClippingDict[read_ID] - 4 >= seq[:30].index(query) + 6:
                                short_tail_outfile.write('{}\t{}\n'.format(read_ID, i))
                                break
                        else: dropped_read2.append([read_ID, 'no_tail'])
                #Standards, not going to be in the soft clipped dict.
                elif read_ID in standards:
                    #All of these reads are in the Read2STAR aligned file.
                    # if read_ID == '1207:14335:91178': pdb.set_trace()
                    firstBasesList = [[k, len(list(g))] for k, g in groupby(seq)]
                    if firstBasesList[0][0] == 'T':
                        TL = firstBasesList[0][1]
                        short_tail_outfile.write('{}\t{}\n'.format(read_ID, TL))
                        num_short_tails += 1
                    elif len(firstBasesList) > 1 and firstBasesList[0][0] == 'A' and firstBasesList[1][0] == 'T': #uridylation
                        TL = firstBasesList[1][1] #Just the number of Ts
                        short_tail_outfile.write('{}\t{}\n'.format(read_ID, TL))
                        num_short_tails += 1
                    elif len(firstBasesList) > 1 and firstBasesList[0][1] <= 2 and firstBasesList[1][0] == 'T': #allow guanylation for fewer than 2 Gs.
                        TL = firstBasesList[1][1] #Just the number of Ts
                        short_tail_outfile.write('{}\t{}\n'.format(read_ID, TL))
                        num_short_tails += 1
                    #look for tails >= 8 and <19 nt anywhere, no mismatches, doesn't need soft clipping as standard. 
                    else:
                        for i in range(19,7,-1):
                            query = ('T'*i)
                            if query in seq[:30]:
                                short_tail_outfile.write('{}\t{}\n'.format(read_ID, i))
                                break
                        else: dropped_read2.append([read_ID, 'no_tail'])
                else:
                    #This statment looks for the original definition of a tail as implemented by SWE in the case that read2 doesn't map
                    #These reads are funneled to the HMM
                    r = regex.compile('(%s){e<=1}' % ('T'*12))
                    match = r.search(seq[:30])
                    if match is not None:
                        match_start = match.start()
                        if seq[match_start] != 'T': match_start += 1 #if error is first pos.
                        if seq[match_start] != 'T': print("Error with starting position.")
                        match_start = match_start + config.TRIM_BASES 
                        new_keep_dict[read_ID] = (keep_dict[read_ID], match_start)
    short_tail_outfile.close()
    return new_keep_dict, dropped_read2, num_short_tails