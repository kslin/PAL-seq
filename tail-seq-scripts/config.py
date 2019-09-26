import numpy as np
import re

'''
converting between Splint and Direct ligation datasets
##Did you make the following changes to the pipeline before running?
1. change the read command for unzipping the input files in the make file
2. change the fastq GZIP to 'else' to read .tar.gz in the config file.
3. change the trim bases to 0 in the config file.
4. source the virtual env.
'''

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

FUTURES = 1 # number of processes, set to 1 if not using multiprocessing
LEN1, LEN2 = 26, 49 # length of read1 and read2
SIGNAL_COL_START = 4 # column in intensity file where intensity values start
SIGNAL_COL_END = (4*(LEN1 + LEN2)) + 4 # column in intensity file where intensity values end
NUM_SKIP, NUM_SKIP_2 = 10, 6 # number of nucleotides to skip in the beginning of read 1 or end of read 2, respectively, when normalizing intensities (usually 10 and 5)
NUC_ORDER = ['A','C','G','T'] # order of nucleotides in intensity file
# UPPERBOUND, LOWERBOUND = 5, -5 # bounds for normalized T-signal
CHUNKSIZE = 1000 # number of intensities to evaluate at a time
NAN_LIMIT = 5 # number of positions with no signal in the intensity file allowed
QUAL = True #If true, only one mismatch is allowed to call a tail with the first 11 nt being T. Else, two mismatches are allowed. 
#Gzip status for fastq and intensity files
FASTQ_GZIP = False #if false, extension is .txt, if True extension is .txt.gz
TRIM_BASES = 14 #adapter length of read 1
OPTIM_CONC = 2 #cycle no corresponding to the conccentration that is used for the linear model. 0 is the background (no streptavidin) cycle

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
