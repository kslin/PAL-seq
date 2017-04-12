import math
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from sklearn.mixture import GaussianMixture

import config


def get_identifier(raw_header):
    return ':'.join(raw_header.split(':')[2:5]).split('#')[0]

def fix_vals(x):
    if x > 0:
        return x
    return 100.0

v_fix_vals = np.vectorize(fix_vals)

def get_normalized_intensities(intensities, read1_sequence, len1, len2):
    """Get average intensities for each nucleotide for use in normalizing"""

    # skip the first 4 elements, which identify the read
    # then skip the number of starting nucleotides specified in the config file
    try:
        intensities = intensities[config.NUM_SKIP:,:]
        read1_sequence = read1_sequence[config.NUM_SKIP:]
    except:
        print("Intensities and sequences for read1 must be longer than config.NUM_SKIP")
        sys.exit()

    assert(len(intensities) == (len1 + len2 - config.NUM_SKIP)), "Intensities length not equal to l1 and l2"
    assert(len(read1_sequence) == len1 - config.NUM_SKIP), "Read1 length must equal l1"

    # convert read1_sequence into one-hot encoding of 4 bits
    read1_sequence = np.array([[int(nt == x) for x in config.NUC_ORDER] for nt in read1_sequence])

    # assert(np.sum(read1_sequence) == (len1 - config.NUM_SKIP)), "Non-{} bases found in read1".format(config.NUC_ORDER)

    # add up the counts for each nucleotide
    read1_nt_counts = np.sum(read1_sequence, axis=0)

    # return None if one or more of the nucleotides doesn't show up in the read1_sequence
    if np.min(read1_nt_counts) == 0:
        return None

    # convert intensities to an array and split into read1 and read2 intensities
    # intensities = np.array([[int(x) for x in i.split()] for i in intensities])
    read1_intensities = intensities[:len1 - config.NUM_SKIP]
    read2_intensities = intensities[-1 * len2:]

    # if any intensities are negative, convert to 1
    read1_intensities = v_fix_vals(read1_intensities)
    read2_intensities = v_fix_vals(read2_intensities)

    # multiply read1 intensities by one-hot to get intensities for the base called
    read1_intensities = np.multiply(read1_sequence, read1_intensities)

    # add up read1_intensities for each nucleotide
    read1_intensities = np.sum(read1_intensities, axis=0)

    # get average read1_intensities for each nucleotide
    norm_vals = np.divide(read1_intensities, read1_nt_counts)

    read2_intensities_normed = np.divide(read2_intensities, norm_vals)

    return read2_intensities_normed

def get_t_signal(intensities, upperbound=config.UPPERBOUND, lowerbound=config.LOWERBOUND):

    # find which index T is at in the intensity file
    t_index = config.NUC_ORDER.index('T')
    other_index = list(range(len(config.NUC_ORDER)))
    other_index.remove(t_index)
    
    # separate T channel from the others
    t_channel = intensities[:, t_index]
    other_channels = intensities[:, other_index]
    background = np.sum(other_channels, axis=1)

    t_signal = np.log2(np.divide(t_channel, background))

    # bound the signal
    t_signal = np.minimum(t_signal, upperbound)
    t_signal = np.maximum(t_signal, lowerbound)

    return t_signal

def get_normalized_intensities2(intensities, read1_sequence, len1, len2):
    """Get average intensities for each nucleotide for use in normalizing"""

    # skip the first 4 elements, which identify the read
    # then skip the number of starting nucleotides specified in the config file
    try:
        intensities = intensities[config.NUM_SKIP:,:]
        read1_sequence = read1_sequence[config.NUM_SKIP:]
    except:
        print("Intensities and sequences for read1 must be longer than config.NUM_SKIP")
        sys.exit()

    assert(len(intensities) == (len1 + len2 - config.NUM_SKIP)), "Intensities length not equal to l1 and l2"
    assert(len(read1_sequence) == len1 - config.NUM_SKIP), "Read1 length must equal l1"

    # convert read1_sequence into one-hot encoding of 4 bits
    read1_sequence = np.array([[int(nt == x) for x in config.NUC_ORDER] for nt in read1_sequence])

    # assert(np.sum(read1_sequence) == (len1 - config.NUM_SKIP)), "Non-{} bases found in read1".format(config.NUC_ORDER)

    # add up the counts for each nucleotide
    read1_nt_counts = np.sum(read1_sequence, axis=0)

    # return None if one or more of the nucleotides doesn't show up in the read1_sequence
    if np.min(read1_nt_counts) == 0:
        return None

    # convert intensities to an array and split into read1 and read2 intensities
    read1_intensities = intensities[:len1 - config.NUM_SKIP]
    read2_intensities = intensities[-1 * len2:]

    # multiply read1 intensities by one-hot to get intensities for the base called
    read1_intensities = np.multiply(read1_sequence, read1_intensities)

    # add up read1_intensities for each nucleotide
    read1_intensities = np.sum(read1_intensities, axis=0)

    # get average read1_intensities for each nucleotide
    norm_vals = np.divide(read1_intensities, read1_nt_counts)

    # intensities_normed = np.divide(intensities, norm_vals)

    return np.divide(intensities, norm_vals)

def get_t_signal2(intensities, len1, len2, upperbound=1.6):

    # find which index T is at in the intensity file
    t_index = config.NUC_ORDER.index('T')
    other_index = list(range(len(config.NUC_ORDER)))
    other_index.remove(t_index)

    t_channel = intensities[:, t_index]
    other_channel = np.max(intensities[:, other_index], axis=1)

    t_signal = t_channel - other_channel

    # fit mixture model to determine normalization constant
    A = GaussianMixture(n_components=2, covariance_type='diag')
    A.fit(np.array(list(np.maximum(t_signal[:2*(len1 - config.NUM_SKIP)],0)) + [0]*10).reshape(-1,1))
    norm_constant = np.max(A.means_)
    if norm_constant <= 0:
        norm_constant = np.max(t_signal[:len1 - config.NUM_SKIP])
    t_signal = np.divide(t_signal, norm_constant)

    t_signal = np.minimum(t_signal, upperbound)

    return t_signal[-1 * len2:]


def get_tail_length2(seq_id, t_signal, a_mu, a_sigma, nota_mu, nota_sigma, rates):
    length = len(t_signal)
    scores = np.zeros((3, length))
    a_scores = norm.logpdf(t_signal, loc=a_mu, scale=a_sigma)
    nota_scores = norm.logpdf(t_signal, loc=nota_mu, scale=nota_sigma)

    prev = 0
    for i in range(length):
        prev += a_scores[i]
        scores[0, i] = prev

    scores[1,0] = -1 * np.inf
    scores[2,0] = nota_scores[0]
    path = np.zeros((2, length), dtype=int)
    rate_dict = {}
    
    
    prev = -1*np.inf
    prev_signal = t_signal[0]
    for j, t in enumerate(t_signal[1:]):
        rate = rates[j+1]
        trans_start = int(max(1, j+2 - max(1, np.round((t - a_mu[j+1])/rate))))
        trans_len = int(j+2 - trans_start)
        ix = np.arange(trans_start, j+2)
        score = np.sum(norm.logpdf(t_signal[ix], loc=(((a_mu[ix] + t_signal[ix])/2) + (rate*(np.arange(1, 1+trans_len)))), scale=a_sigma))
        continue_transition_score = scores[0, trans_start - 1] + score
        new_transition_score = scores[0,j] + norm.logpdf(t, loc=a_mu[j+1]+rate, scale=a_sigma)

        if (continue_transition_score > new_transition_score):
            new_rate = (t - t_signal[trans_start-1])/trans_len
            if new_rate < 0:
                rate_dict[j+1] = new_rate
            else:
                rate_dict[j+1] = np.nan
            path[0, j+1] = trans_len
            new_score = continue_transition_score
        else:
            new_rate = t - t_signal[j]
            if new_rate < 0:
                rate_dict[j+1] = new_rate
            else:
                rate_dict[j+1] = np.nan
            path[0, j+1] = 1
            new_score = new_transition_score
            
        if prev > scores[2,j]: # T -> N
            scores[2, j+1] = prev + nota_scores[j+1]
            path[1, j+1] = 1
        else:
            scores[2, j+1] = scores[2,j] + nota_scores[j+1]

        scores[1, j+1] = new_score
        prev = new_score
        prev_signal = t
    
    start = np.argmax(scores[:, -1])
    if start == 0:
        return seq_id, length, t_signal, [np.nan]*length, [np.nan]*length
    elif start == 1:
        trans_length = path[0,-1]
        tail_length = int(length - trans_length)
        return seq_id, tail_length, list(t_signal[:tail_length]) + [np.nan]*trans_length, [np.nan]*length, [np.nan]*(tail_length) + [rate_dict[length-1]]*trans_length
    else:
        if 1 not in list(path[1,:]):
            return seq_id, 0, [np.nan]*length, t_signal, [np.nan]*length
        n_length = list(path[1,:])[::-1].index(1)
        trans_length = int(path[0,(length - n_length - 1)])
        tail_length = int(length - n_length - trans_length)
        A_vals = list(t_signal[:tail_length]) + [np.nan]*(trans_length + n_length)
        notA_vals = [np.nan]*(tail_length + trans_length) + list(t_signal[tail_length + trans_length:])
        new_rates = [np.nan]*(tail_length) + [rate_dict[(length - n_length - 1)]]*trans_length + [np.nan]*n_length
        return seq_id, tail_length, A_vals, notA_vals, new_rates

def get_tail_length3(t_signal, a_mu, a_sigma, nota_mu, nota_sigma, rates):
    length = len(t_signal)
    a_scores = norm.logpdf(t_signal, loc=a_mu, scale=a_sigma)
    nota_scores = norm.logpdf(t_signal, loc=nota_mu, scale=nota_sigma)

    scores = np.zeros((3, length))
    prev = 0
    for i in range(length):
        prev += a_scores[i]
        scores[0, i] = prev

    scores[1,0] = -1 * np.inf
    scores[2,0] = nota_scores[0]

    path = np.zeros((2, length), dtype=int)
    
    prev_score = -1*np.inf
    prev_num = 0
    prev_signal = t_signal[0]
    for j, t in enumerate(t_signal[1:]):
        rate = rates[j+1]
        # prob of transition starting here
        new_transition_score = scores[0, j] + norm.logpdf(t, loc=a_mu[j+1]+rate, scale=a_sigma)

        # prob of continuing transition state
        continue_transition_score = prev_score + norm.logpdf(t, loc=a_mu[j+1]+((prev_num+1)*rate), scale=a_sigma)
        print(j+1)
        print(scores[0, j], scores[0, j+1])
        print(scores[1, j], new_transition_score, continue_transition_score)
        print(scores[2, j], nota_scores[j+1])

        if (continue_transition_score < new_transition_score):
            path[0, j+1] = 1
            new_score = new_transition_score
            prev_num = 1
        else:
            path[0, j+1] = prev_num + 1
            new_score = continue_transition_score
            prev_num += 1
            
        if prev_score > scores[2,j]: # T -> N
            scores[2, j+1] = prev_score + nota_scores[j+1]
            path[1, j+1] = 1
        else: # N -> N
            scores[2, j+1] = scores[2,j] + nota_scores[j+1]

        scores[1, j+1] = new_score
        prev_score = new_score
        prev_signal = t
    
    start = np.argmax(scores[:, -1])
    if start == 0:
        tail_length = length
        A_vals = t_signal
        notA_vals = [np.nan]*length
        rate_vals = [np.nan]*length
    elif start == 1:
        trans_length = path[0,-1]
        tail_length = int(length - trans_length)
        A_vals = list(t_signal[:tail_length]) + [np.nan]*trans_length
        notA_vals = [np.nan]*length
        rate_vals = [np.nan]*(tail_length) + list((t_signal[tail_length:] - a_mu[tail_length:]) / np.arange(1,trans_length + 1))
    else:
        if 1 not in list(path[1,:]):
            tail_length = 0
            A_vals = [np.nan]*length
            notA_vals = t_signal
            rate_vals = [np.nan]*length
        else:
            n_length = list(path[1,:])[::-1].index(1) + 1
            trans_length = int(path[0, (length - n_length - 1)])
            tail_length = int(length - n_length - trans_length)
            trans_index, n_index = tail_length, tail_length + trans_length

            A_vals = list(t_signal[:tail_length]) + [np.nan]*(trans_length + n_length)
            notA_vals = [np.nan]*(n_index) + list(t_signal[n_index:])
            rate_vals = [np.nan]*(tail_length) + \
                        list((t_signal[trans_index:n_index] - a_mu[trans_index:n_index]) / np.arange(1,trans_length + 1)) + \
                        [np.nan]*n_length
    
    return tail_length, np.array(A_vals).reshape((1,length)), np.array(notA_vals).reshape((1,length)), np.array(rate_vals).reshape((1,length))


def get_tail_length_batch(seq_ids, t_signals, a_mu, a_sigma, nota_mu, nota_sigma, rates):
    length = len(a_mu)

    tail_lengths, a_val_arrays, nota_val_arrays, rate_arrays = [], [], [], []
    for i in range(t_signals.shape[0]):
        t_signal = t_signals[i, :]
        tail_length, A_vals, notA_vals, rate_vals = get_tail_length3(t_signal, a_mu, a_sigma, nota_mu, nota_sigma, rates)
        tail_lengths.append(tail_length)
        a_val_arrays.append(A_vals)
        nota_val_arrays.append(notA_vals)
        rate_arrays.append(rate_vals)

    return list(seq_ids), tail_lengths, a_val_arrays, nota_val_arrays, rate_arrays


def add_logpdfs(x, y, min_val=-700):
    if (x < min_val) & (y < min_val):
        return min_val
    if abs(x-y) > 4:
        return max(x,y)
    return min(x,y) + np.log1p(np.exp(abs(x-y)))

add_logpdfs_vec = np.vectorize(add_logpdfs)


def exp_fit(x, k):
    return np.exp(((-1*k )/ (x))) - 1


def get_tail_length_gmm(t_signal, tail_start, a_mu, a_sigma, nota_mu1, nota_sigma1, nota_mu2, nota_sigma2, ratio, rates):

    original_length = len(t_signal)
    t_signal = t_signal[tail_start:]
    a_mu = a_mu[tail_start:]
    a_sigma = a_sigma[tail_start:]
    nota_mu1 = nota_mu1[tail_start:]
    nota_sigma1 = nota_sigma1[tail_start:]
    rates = rates[tail_start:]

    length = len(t_signal)
    a_scores = norm.logpdf(t_signal, loc=a_mu, scale=a_sigma)
    # nota_scores1 = norm.logpdf(t_signal, loc=nota_mu1, scale=nota_sigma1) + np.log(ratio)
    # nota_scores2 = norm.logpdf(t_signal, loc=nota_mu2, scale=nota_sigma2) + np.log(1.0 - ratio)

    # nota_scores = add_logpdfs_vec(nota_scores1,nota_scores2)

    nota_scores = norm.logpdf(t_signal, loc=nota_mu1, scale=nota_sigma1)

    scores = np.zeros((3, length))
    prev = 0
    for i in range(length):
        prev += a_scores[i]
        scores[0, i] = prev

    scores[1,0] = -1 * np.inf
    scores[2,0] = nota_scores[0]

    path = np.zeros((2, length), dtype=int)
    
    prev_score = -1*np.inf
    prev_num = 0
    prev_signal = t_signal[0]
    for j, t in enumerate(t_signal[1:]):
        rate = rates[j+1]
        # prob of transition starting here
        new_transition_score = scores[0, j] + norm.logpdf(t, loc=a_mu[j+1]+rate, scale=a_sigma[j+1])

        # prob of continuing transition state
        continue_transition_score = prev_score + norm.logpdf(t, loc=a_mu[j+1]+((prev_num+1)*rate), scale=a_sigma[j+1])

        if (continue_transition_score < new_transition_score):
            path[0, j+1] = 1
            new_score = new_transition_score
            prev_num = 1
        else:
            path[0, j+1] = prev_num + 1
            new_score = continue_transition_score
            prev_num += 1
            
        if prev_score > scores[2,j]: # T -> N
            scores[2, j+1] = prev_score + nota_scores[j+1]
            path[1, j+1] = 1
        else: # N -> N
            scores[2, j+1] = scores[2,j] + nota_scores[j+1]

        scores[1, j+1] = new_score
        prev_score = new_score
        prev_signal = t
    
    start = np.argmax(scores[:, -1])
    if start == 0:
        tail_length = length
        A_vals = list(t_signal)
        notA_vals = [np.nan]*length
        rate_vals = [np.nan]*length
    elif start == 1:
        trans_length = path[0,-1]
        tail_length = int(length - trans_length)
        A_vals = list(t_signal[:tail_length]) + [np.nan]*trans_length
        notA_vals = [np.nan]*length
        rate_vals = [np.nan]*(tail_length) + list((t_signal[tail_length:] - a_mu[tail_length:]) / np.arange(1,trans_length + 1))
    else:
        if 1 not in list(path[1,:]):
            tail_length = 0
            A_vals = [np.nan]*length
            notA_vals = list(t_signal)
            rate_vals = [np.nan]*length
        else:
            n_length = list(path[1,:])[::-1].index(1) + 1
            trans_length = int(path[0, (length - n_length - 1)])
            tail_length = int(length - n_length - trans_length)
            trans_index, n_index = tail_length, tail_length + trans_length

            A_vals = list(t_signal[:tail_length]) + [np.nan]*(trans_length + n_length)
            notA_vals = [np.nan]*(n_index) + list(t_signal[n_index:])
            rate_vals = [np.nan]*(tail_length) + \
                        list((t_signal[trans_index:n_index] - a_mu[trans_index:n_index]) / np.arange(1,trans_length + 1)) + \
                        [np.nan]*n_length
    
    return tail_length, np.array(([np.nan]*tail_start) + A_vals).reshape((1,original_length)), np.array(([np.nan]*tail_start) + notA_vals).reshape((1,original_length)), np.array(([np.nan]*tail_start) + rate_vals).reshape((1,original_length))


def get_tail_length_batch_gmm(seq_ids, starts, t_signals, a_mu, a_sigma, nota_mu1, nota_sigma1, nota_mu2, nota_sigma2, ratio, rates):
    length = len(a_mu)

    tail_lengths, a_val_arrays, nota_val_arrays, rate_arrays = [], [], [], []
    for i in range(t_signals.shape[0]):
        t_signal = t_signals[i, :]
        tail_length, A_vals, notA_vals, rate_vals = get_tail_length_gmm(t_signal, starts[i], a_mu, a_sigma, nota_mu1, nota_sigma1,
                                                                        nota_mu2, nota_sigma2, ratio, rates)
        tail_lengths.append(tail_length)
        a_val_arrays.append(A_vals)
        nota_val_arrays.append(notA_vals)
        rate_arrays.append(rate_vals)

    return list(seq_ids), tail_lengths, a_val_arrays, nota_val_arrays, rate_arrays


# a_mu = np.array([1,1,1,1,2,2,2,2,2,2,2,2,2])
# a_sigma = 0.1
# nota_mu1 = a_mu * 0.5
# nota_sigma1 = 0.1
# nota_mu2 = 0
# nota_sigma2 = 0.1
# # nota_sigma1 = np.array([0.1]*10)
# rates = np.array([-1,-1,-1,-0.9,-0.7,-0.6,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5])
# ratio = 0.5

# t_signal = np.array([[1,1,1,0,1,1,0,1,0,1,0,1,1],[1,1,1,1,2,1.5,1,0.5,0,1,1,0,1]])
# t0 = time.time()
# for i in range(t_signal.shape[0]):
#     print(t_signal[i, :])
#     print(get_tail_length_gmm(t_signal[i,:], 0, a_mu, a_sigma, nota_mu1, nota_sigma1, nota_mu2, nota_sigma2, ratio, rates))
# print(time.time() - t0)




