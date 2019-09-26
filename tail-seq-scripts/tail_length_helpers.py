import gzip
import numpy as np
import pandas as pd #needs 0.24 or higher
import config
import pdb
import scipy.stats as ss
import sys

def read_training_set(infile):
    """
    Get standard data for linear model.
    Return a dictionary of numpy arrays, one entry for each standard class.
    The numpy arrays have the first column as the length value and then the different concentrations
    """

    signal_file = pd.read_csv(infile, usecols = [1,*range(3,3+config.NUM_SKIP_2)],sep='\t', header=None, engine='c')

    stds_uniq = signal_file[1].unique()
    training_array_dict = {}
    for std in stds_uniq:
        training_array_dict[std] = signal_file[signal_file[1] == std].to_numpy()

    return training_array_dict


def train_model(training_array_dict, training_param_file, med_param_file):
    """
    Get linear model parameters.
    """
    
    medList = [] #median dict
    for std, arr in training_array_dict.items():
        medList.append(np.median(arr,axis = 0))
    medArr = np.vstack(medList)

    #linear regression
    LinRegTup = [ss.linregress(medArr[:,(0,(1+i))])[:] for i in range(config.NUM_SKIP_2)]
    #slope, intercept, r_value, p_value, std_err
    LinRegArr = np.array(LinRegTup)
    np.savetxt(med_param_file, medArr, header = "\t".join(['nucleotide_length', 'intensity_conc1', 'intensity_conc2', 'intensity_conc3', 'intensity_conc4','intensity_conc5']),comments = '',delimiter = '\t') 
    np.savetxt(training_param_file, LinRegArr, header = "\t".join(['slope', 'intercept', 'r_value', 'p_value', 'std_err']),comments = '',delimiter = '\t') 
    sys.exit()