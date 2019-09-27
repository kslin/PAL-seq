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

    signal_file = pd.read_csv(infile, usecols = [0,1,*range(3,3+config.NUM_SKIP_2)],sep='\t', header=None, engine='c')

    stds_uniq = signal_file[1].unique()
    training_array_dict = {}
    for std in stds_uniq:
        training_array_dict[std] = signal_file[signal_file[1] == std].loc[:,1:].to_numpy()

    return training_array_dict, signal_file.loc[:,(0,1)]


def train_model(training_array_dict, training_param_file, med_param_file):
    """
    Get linear model parameters.
    """
    
    medList = [] #median dict
    for std, arr in training_array_dict.items():
        medList.append(np.median(arr,axis = 0))
    medArr = np.vstack(medList)

    #Changes the tail lengths to account for the cycles prior to extension.
    medArr[:,0] = medArr[:,0] - (config.LEN1 - config.TRIM_BASES)

    #removes negative tail lengths (i.e. 10mer)
    # mask = np.ones(medArr.shape[0], dtype=bool)
    # mask[np.where(medArr[:,0] < 0)] = False
    # medArr = medArr[mask,:]

    #removes the 324mer, as per Alex's notes
    mask = np.ones(medArr.shape[0], dtype=bool)
    mask[np.where(medArr[:,0] == (324 - (config.LEN1 - config.TRIM_BASES)))] = False
    medArr = medArr[mask,:]

    #linear regression
    LinRegTup = [ss.linregress(medArr[:,(0,(1+i))])[:] for i in range(config.NUM_SKIP_2)]

    #slope, intercept, r_value, p_value, std_err
    LinRegArr = np.array(LinRegTup)
    np.savetxt(med_param_file, medArr, header = "\t".join(['nucleotide_length'] + ['conc_' + str(i) for i in range(config.NUM_SKIP_2)]),comments = '',delimiter = '\t') 
    np.savetxt(training_param_file, LinRegArr, header = "\t".join(['slope', 'intercept', 'r_value', 'p_value', 'std_err']),comments = '',delimiter = '\t') 

    return(LinRegArr)

def get_batch_tail_length(signal_file, LinRegArr):
    """
    Given a signal file, use the model to calculate tail lengths
    """
    # extract metadata
    ids, lengths = [], []

    with open(signal_file, 'r') as sf:
        for line in sf:
            read_ID = line.split("\t")[0]
            preExtTl = int(line.split("\t")[1]) #nucleotides of tail prior to extension
            signal = float(line.split("\t")[2 + config.OPTIM_CONC])

            #calculate tail length
            tl = (signal - LinRegArr[config.OPTIM_CONC,1]) / LinRegArr[config.OPTIM_CONC,0] + preExtTl

            ids.append(read_ID)
            lengths.append(tl)
    
    return ids, lengths