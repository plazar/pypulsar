"""
mydetrend.py

A detrending algorithm similar to scipy.signal.detrend. My version
accepts a mask, which can be used to omit values when detrending.

Patrick Lazarus, Nov. 13, 2010
"""

import numpy as np
import scipy.linalg

import matplotlib.pyplot as plt


def detrend(ydata, xdata=None, mask=None, order=1):
    """Detrend 'data' using a polynomial of given order.
    
        Inputs:
            ydata: A 1D array to be detrended.
            xdata: A 1D array of x-values to use
                    (Optional. Default: Use indices at xdata).
            mask: A 1D array of booleans to used to omit
                    values in data when detrending.
                    True means omit value.
                    (Default: Mask nothing)
            order: Order of polynomial to use (Default: 1)

        Output:
            detrended: a 1D array.
    """
    if xdata is None:
        xdata = np.arange(ydata.size)
    powers = np.arange(order+1)

    A = np.repeat(xdata, order+1)
    A.shape = (xdata.size, order+1)
    A = A**powers

    if mask is None:
        unmasked = np.ones(ydata.size, dtype='bool')
    else:
        unmasked = np.bitwise_not(mask) 
    x, resids, rank, s = scipy.linalg.lstsq(A[unmasked], ydata[unmasked])
    detrended =  ydata - np.dot(A, x)
    return detrended

