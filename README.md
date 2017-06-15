# Tail-Seq

Requirements:
- Python 2.7 or above
- ghmm
- pandas
- numpy
- concurrent (part of the standard library in Python 3)

1. Calculate the normalized T-signal from the raw intensity file:

python tail-seq-scripts/get_signal.py [--f1 fastq1.txt.gz] [--l1 read1_length] [--f2 fastq2.txt.gz] [--l2 read2_length] [-i intensityfile.txt.gz] [-a read1_annotations.txt] [-o output_directory] [-f (optional, use -f to multi-thread)]

2. Train Gaussian HMM and call tail-lengths:

python tail-seq-scripts/tail_length_hmm.py [-s normalized_signal.txt (from get_signal.py)] [-l read2_length] [-o output_directory] [-t train_size (optional, default 10000)] [-m max_iter (optional, default 10000)] [--tol tolerance (optional, default 0.01)] 