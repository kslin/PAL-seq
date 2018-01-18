# Tail-Seq

Requirements:
- Python 2.7 or above
- [GHMM](http://ghmm.org/)
- concurrent (part of the standard library in Python 3)
- numpy
- pandas
- matplotlib


1. Calculate the normalized T signal from the raw intensity file:

python tail-seq-scripts/get_signal.py [--f1 fastq1.txt.gz] [--l1 read1 length] [--f2 fastq2.txt.gz] [--l2 read2 length] [-i intensityfile.txt.gz] [-a read1_annotations.txt] [-o output directory] [-f specify number of threads, default 1]

2. Plot T signals

python tail-seq-scripts/plot_t_signals.py [-s normalized_signal.txt (from get_signal.py)] [-o output_directory] [-b number of bins] [-l length of read2]

3. Train Gaussian Mixture HMM and call tail-lengths:

python tail-seq-scripts/tail_length_hmm.py [-s normalized_signal.txt (from get_signal.py)] [-l read2 length] [-o output_directory] [--twostate toggles between 2-state and 3-state model, default False] [-t train size, default 10000] [-m max iter, default 10000)] [--tol tolerance, default 0.01] [-f specify number of threads, default 1]