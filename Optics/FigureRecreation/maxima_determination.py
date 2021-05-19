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
def get_extrema(trace, bias):
    tr_peaks_only = list(trace[find_peaks(trace)[0]])
    tr_peaks_only += list(trace[find_peaks(-trace)[0]])
    tr_mean = np.mean(np.convolve(trace, [0.2, 0.2, 0.2, 0.2, 0.2], mode='same'))
    if len(tr_peaks_only) >= 1:
        groups = group_values(tr_peaks_only, bias)
        vals = [np.mean(gr)-tr_mean for gr in groups]
        return vals
    return [0.0]
