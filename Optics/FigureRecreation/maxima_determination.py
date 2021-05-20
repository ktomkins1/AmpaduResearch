import numpy as np
from scipy.signal import find_peaks

'''
    Group a list based on the distance between values

    param y: the list of values to bin
    param bias: the maximum change between values in each group

    returns: list of groups
'''
def group_values(y, bias):
    grs = []
    prk = -np.inf
    for k in sorted(y):
        if k - bias > prk: grs.append([])
        prk = k
        grs[-1].append(k)
    return grs

'''
    Get the extrema of a (periodic) trace

    param trace: the array of values to calculate on
    param bias: the maximum allowable change in value

    returns: list of extrema, averaged to reduce numerocity
'''
def get_extrema(trace, config):
    bias = config['ex_bias']
    norm = config['bf_norm']
    tr_peaks_only = list(trace[find_peaks(trace)[0]])
    tr_peaks_only += list(trace[find_peaks(-trace)[0]])
    tr_mean = np.mean(np.convolve(trace, [0.2, 0.2, 0.2, 0.2, 0.2], mode='same'))
    if 'bf_norm_type' not in config.keys():
        try:
            tr_mean = config['bf_norm_val']
        except KeyError:
            tr_mean = 1.0
            config['bf_norm_val'] = tr_mean
    if len(tr_peaks_only) >= 1:
        groups = group_values(tr_peaks_only, bias)
        vals = []
        for gr in groups:
            pt = np.mean(gr)
            if norm: pt -= tr_mean
            vals.append(pt)
        return vals
    if norm:
        return [np.mean(trace)-tr_mean]
    else:
        return [tr_mean]
