import numpy as np

def fastq_header_to_ID(line):
  """
  Convert fastq header to read ID
  """
  return line.split('#')[0].replace('WIGTC-HISEQ:2:','')

def intensity_line_to_ID(line):
  """
  Extract read ID from intensity line
  """
  read_ID = ':'.join(line[:3])
  # read_ID = '{}:{}:{}#{}'.format(line[0], line[1], line[2], line[3])
  return read_ID

### data configurations ###

LEN1, LEN2 = 50, 250

SIGNAL_COL_START = 4
SIGNAL_COL_END = (4*(LEN1 + LEN2)) + 4

NUM_SKIP = 10

NUC_ORDER = ['A','C','G','T']

UPPERBOUND, LOWERBOUND = 5, -5

CHUNKSIZE = 1000

NAN_LIMIT = 5

WINDOW_SIZE, CUTOFF = 40, -20

START_SIGNAL = 100.0

TRAIN_SIZE = 10000
MAX_ITER = 10000
TOL = 0.01

HMM_PARAMS_3STATE = {

                        "TRANSITION_INIT": np.array([[0.95, 0.03, 0.019, 0.001],
                                                     [0.001, 0.749, 0.249, 0.001],
                                                     [0.001, 0.001, 0.997, 0.001],
                                                     [0.95, 0.025, 0.0249, 0.001]]),

                        "MEANS_INIT": np.array([[1.0], [0.0], [-1.0], [START_SIGNAL]]),

                        "VARS_INIT": np.array([[0.25], [0.25], [1.0], [1.0]]),

                        "START_PROB_INIT": np.array([0.001, 0.001, 0.001, 0.997])
                      }


HMM_PARAMS_2STATE = {
                        "TRANSITION_INIT": np.array([[0.95, 0.049, 0.001],
                                                    [0.001, 0.998, 0.001],
                                                    [0.95, 0.049, 0.001]]),

                        "MEANS_INIT": np.array([[1.0], [-1.0], [START_SIGNAL]]),

                        "VARS_INIT": np.array([[0.25], [0.25], [1.0]]),

                        "START_PROB_INIT": np.array([0.001, 0.001, 0.998])
                      }

# HMM_PARAMS_2STATE = {
#                         "TRANSITION_INIT": np.array([[0.95, 0.05, 0.0],
#                                                     [0.0, 1.0, 0.0],
#                                                     [0.95, 0.05, 0.0]]),

#                         "MEANS_INIT": np.array([[1.0], [-1.0], [START_SIGNAL]]),

#                         "VARS_INIT": np.array([[0.25], [0.25], [1.0]]),

#                         "START_PROB_INIT": np.array([0.0, 0.0, 1.0])
#                       }
