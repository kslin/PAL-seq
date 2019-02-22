import numpy as np
import re

### Functions for extracting the same read ID's from fastq files and intensity files ###

def fastq_header_to_ID(line):
  """
  Convert fastq header to read ID
  """
  searchSTR = re.compile('WIGTC-HISEQ:\\d:')
  return re.sub(searchSTR,'',line.split('#')[0])

def intensity_line_to_ID(line):
  """
  Extract read ID from intensity line
  """
  read_ID = ':'.join(line[:3])
  return read_ID

### Preprocess run configurations ###

FUTURES = 24 # number of processes, set to 1 if not using multiprocessing
LEN1, LEN2 = 50, 250 # length of read1 and read2
SIGNAL_COL_START = 4 # column in intensity file where intensity values start
SIGNAL_COL_END = (4*(LEN1 + LEN2)) + 4 # column in intensity file where intensity values end
NUM_SKIP = 10 # number of nucleotides to skip when normalizing intensities
NUC_ORDER = ['A','C','G','T'] # order of nucleotides in intensity file
UPPERBOUND, LOWERBOUND = 5, -5 # bounds for normalized T-signal
CHUNKSIZE = 1000 # number of intensities to evaluate at a time
NAN_LIMIT = 5 # number of positions with no signal in the intensity file allowed
QUAL = True #If true, only one mismatch is allowed to call a tail with the first 11 nt being T. Else, two mismatches are allowed. 
#Gzip status for fastq and intensity files
FASTQ_GZIP = "None" #if false, extension is .txt, if False extension is .txt.gz, else extension is .tar.gz
TRIM_BASES = 4 #Should be set to 0 for a splint ligation, 4 or 8 for a direct lig run. 
# INTENSITY_GZIP = False #if false, extension is .txt, else extension is .txt.gz ##This is handled differently in the get_signal_helpers code


# WINDOW_SIZE, CUTOFF = 40, -20

### HMM configurations ###

START_SIGNAL = 100.0 # HMM value for start of tail
TRAIN_SIZE = 10000 # number of reads to train HMM on
MAX_ITER = 10000 # maximum iterations for training HMM
TOL = 0.01 # tolerance for training HMM

# initial parameters for 3-state HMM
HMM_PARAMS_3STATE = {

                        "TRANSITION_INIT": np.array([[0.95, 0.03, 0.019, 0.001],
                                                     [0.001, 0.749, 0.249, 0.001],
                                                     [0.001, 0.001, 0.997, 0.001],
                                                     [0.95, 0.025, 0.024, 0.001]]),

                        "MEANS_INIT": np.array([[1.0], [0.0], [-1.0], [START_SIGNAL]]),

                        "VARS_INIT": np.array([[0.25], [0.25], [1.0], [1.0]]),

                        "START_PROB_INIT": np.array([0.001, 0.001, 0.001, 0.997])
                      }

# initial parameters for 2-state HMM
HMM_PARAMS_2STATE = {
                        "TRANSITION_INIT": np.array([[0.95, 0.049, 0.001],
                                                    [0.001, 0.998, 0.001],
                                                    [0.95, 0.049, 0.001]]),

                        "MEANS_INIT": np.array([[1.0], [-1.0], [START_SIGNAL]]),

                        "VARS_INIT": np.array([[0.25], [0.25], [1.0]]),

                        "START_PROB_INIT": np.array([0.001, 0.001, 0.998])
                      }
