import gzip
import numpy as np
import pandas as pd #needs 0.24 or higher
import config
import pdb
import scipy.stats as ss
import sys
import statsmodels.nonparametric.kernel_regression as kr
import nlopt

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

def optimize_logistic(params, grad, x, data):
    params = np.exp(params) #exponentiate
    y = gen_logistic(params, x)
    SUM_OF_SQUARES = np.sum((data - y)**2)
    return(SUM_OF_SQUARES)

def gen_logistic(params, x):
    A = params[0] #lower asymptote
    L = params[1] #max val
    k = params[2] #steepness
    m = params[3] #x val at midpoint

    y = L / (1 + np.exp(-k * (x - m))) - A

    return(y)

def gen_logistic_rev(params, y):
    A = params[0] #lower asymptote
    L = params[1] #max val
    k = params[2] #steepness
    m = params[3] #x val at midpoint

    x = np.log((L - (y + A)) / (y + A)) / -k + m

    return(x)

def train_model(training_array_dict, training_param_file, med_param_file):
    """Get linear model parameters.
    (Or LOESS interpolation.)
    Keyword arguments
        training_array_dict: dictionary of intensities for each standard
        training_param_file: output file path for the training parameters
        med_param_file: output file path for the normalized intensities
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

    #nlopt logistic
    opt = nlopt.opt(nlopt.LN_NELDERMEAD, 4)
    opt.set_min_objective(lambda params, grad: optimize_logistic(params, grad, x = medArr[:,0][:], data = medArr[:,config.OPTIM_CONC + 1][:]))
    opt.set_ftol_rel(1E-5)
    xopt = np.exp(opt.optimize(np.log(np.array([5, 1, 0.1, 100]))))
    opt_val = opt.last_optimum_value()
    result = opt.last_optimize_result()

    #linear regression
    #Note on 2020 04 09: I've made the dependent variable the std tail length
    # LinRegTup = [ss.linregress(medArr[:,((0, 1+i))])[:] for i in range(config.NUM_SKIP_2)]
    LinRegTup = []
    for i in range(config.NUM_SKIP_2):
        LinRegTup.append(ss.linregress(medArr[:,((1+i), 0)])[:])
    #Kernel regression.
    # KrObj = kr.KernelReg(medArr[:,0:1], medArr[:,2:3], 'c', 'll')

    #slope, intercept, r_value, p_value, std_err
    LinRegArr = np.array(LinRegTup)
    np.savetxt(med_param_file, medArr, header = "\t".join(['nucleotide_length'] + ['conc_' + str(i) for i in range(config.NUM_SKIP_2)]),comments = '',delimiter = '\t') 
    np.savetxt(training_param_file, LinRegArr, header = "\t".join(['slope', 'intercept', 'r_value', 'p_value', 'std_err']),comments = '',delimiter = '\t') 

    #Return only the values from the optim_conc data
    return(LinRegArr[config.OPTIM_CONC], xopt)

def get_batch_tail_length(signal_file, LinRegArr):
    """
    Given a signal file, use the linear model to calculate tail lengths
    """
    # extract metadata
    ids, lengths, signal_arr, preExtT_arr = [], [], [], []
    slope  = LinRegArr[0]
    intercept = LinRegArr[1]

    with open(signal_file, 'r') as sf:
        for index, line in enumerate(sf):
            read_ID = line.split("\t")[0]
            preExtTl = int(line.split("\t")[1]) #nucleotides of tail prior to extension
            signal = float(line.split("\t")[2 + config.OPTIM_CONC])
            #calculate tail length
            #original
            # tl = (signal - LinRegArr[config.OPTIM_CONC,1]) / LinRegArr[config.OPTIM_CONC,0]
            tl = signal * slope + intercept
            tl = np.max((0, tl)) + preExtTl #ReLU tl values. 
            lengths.append(tl)
            ids.append(read_ID)

        ## These lines run the kernel regression method. 
        #     signal_arr.append(signal)
        #     preExtT_arr.append(preExtTl)
        #     if index != 0 and index % 10000 == 0:
        #         tl = KrObj.fit(np.array(signal_arr))[0] + np.array(preExtT_arr)
        #         lengths += tl.tolist()
        #         signal_arr = []
        #         preExtT_arr = []
        # tl = KrObj.fit(np.array(signal_arr))[0] + np.array(preExtT_arr)
        # lengths += tl.tolist()

    return ids, lengths

def get_batch_tail_length_logistic(signal_file, LogitParams):
    """
    Given a signal file, use the logistic model to calculate tail lengths
    """
    # extract metadata
    ids, lengths, signal_arr, preExtT_arr = [], [], [], []

    with open(signal_file, 'r') as sf:
        for index, line in enumerate(sf):
            read_ID = line.split("\t")[0]
            preExtTl = int(line.split("\t")[1]) #nucleotides of tail prior to extension
            signal = float(line.split("\t")[2 + config.OPTIM_CONC])
            #calculate tail length
            # pdb.set_trace()
            tl = gen_logistic_rev(LogitParams, signal)
            tl = np.nan_to_num(tl) + preExtTl #remove nan
            lengths.append(tl)
            ids.append(read_ID)

    return ids, lengths