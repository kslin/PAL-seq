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
  read_ID = '{}:{}:{}'.format(line[0], line[1], line[2])
  return read_ID

### data configurations ###

LEN1, LEN2 = 50, 250

SIGNAL_COL_START = 4
SIGNAL_COL_END = (4*(LEN1 + LEN2)) + 4

NUM_SKIP = 10

NUC_ORDER = ['A','C','G','T']

UPPERBOUND, LOWERBOUND = 5, -5

CHUNKSIZE = 5000

NAN_LIMIT = 5

WINDOW_SIZE, CUTOFF = 40, -20

START_SIGNAL = 10000000000.0

HMM_PARAMS_3STATE = {
                        "TRANSITIONMATRIX": [[0.94, 0.03, 0.01, 0.01, 0.01],
                                             [0.0, 0.5, 0.4, 0.08, 0.02],
                                             [0.0, 0.0, 0.6, 0.38, 0.02],
                                             [0.0, 0.0, 0.0, 0.95, 0.05],
                                             [0.95, 0.01, 0.01, 0.03, 0.0]],

                        "EMISSIONMATRIX": [[[1.5, 0.0], [1.5, 1.5], [0.95, 0.05]],
                                           [[1.5, -1.0], [1.5, 1.5], [0.75, 0.25]],
                                           [[1.5, -1.0], [1.5, 1.5], [0.5, 0.5]],
                                           [[1.5, -1.0], [1.5, 1.5], [0.25, 0.75]],
                                           [[START_SIGNAL, 0.0], [1.0, 1.0], [1.0, 0.0]]],

                        "PI": [0.01, 0.01, 0.01, 0.01, 0.96]
                      }


HMM_PARAMS_2STATE = {
                      "TRANSITIONMATRIX": [[0.95,0.05,0.0],
                                           [0.0,1.0,0.0],
                                           [0.95,0.05,0.0]],

                      "EMISSIONMATRIX": [[1.0,0.5],
                                         [-1.0,0.5],
                                         [10000000000.0,1.0]],

                      "PI": [0.0, 0.0, 1.0]
                      }
