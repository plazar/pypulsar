"""
mydetrend.py

A detrending algorithm similar to scipy.signal.detrend. My version
accepts a mask, which can be used to omit values when detrending.

Patrick Lazarus, Nov. 13, 2010
"""

import numpy as np
import scipy.linalg

def detrend(data, order=1, mask=None):
    """Detrend 'data' using a polynomial of given order.
    
        Inputs:
            data: a 1D array to be detrended.
            order: Order of polynomial to use (Default: 1)
            mask: A 1D array of booleans to used to omit
                    values in data when detrending.
                    True means omit value.
                    (Default: Mask nothing)

        Output:
            detrended: a 1D array.
    """
    indices = np.arange(data.size)
    powers = np.arange(order+1)

    A = np.repeat(np.arange(indices.size), order+1)
    A.shape = (indices.size, order+1)
    A = A**powers

    if mask is None:
        unmasked = np.ones(data.size, dtype='bool')
    else:
        unmasked = np.bitwise_not(mask) 
    x, resids, rank, s = scipy.linalg.lstsq(A[unmasked], data[unmasked])
    detrended =  data - np.dot(A, x)
    return detrended


