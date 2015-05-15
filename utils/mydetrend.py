"""
mydetrend.py

A detrending algorithm similar to scipy.signal.detrend. My version
accepts a mask, which can be used to omit values when detrending.

Patrick Lazarus, Nov. 13, 2010
"""
import datetime

import numpy as np
import scipy.linalg

import matplotlib.pyplot as plt

DEBUG = 0

def old_detrend(ydata, xdata=None, mask=None, order=1):
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
    print x
    detrended =  ydata - np.dot(A, x)
    if DEBUG:
        fn = "mydetrend_%s.txt" % datetime.datetime.isoformat(datetime.datetime.now())
        f = open(fn, 'w')
        for val in ydata:
            f.write("%.15g\n" % val)
        f.close()
        plt.plot(ydata, label="Input data")
        plt.plot(np.dot(A, x), label="Order %d polynomial" % order)
        plt.plot(detrended, label="Detrended data")
        plt.legend(loc='best', prop=dict(size='x-small'))
        plt.show()
    return detrended


def detrend(ydata, xdata=None, order=1, bp=[], numpieces=None):
    """Detrend 'data' using a polynomial of given order.
    
        Inputs:
            ydata: A 1D array to be detrended.
            xdata: A 1D array of x-values to use
                (Default: Use indices at xdata).
            order: Order of polynomial to use (Default: 1)
            bp: Breakpoints. Break the input data into segments
                that are detrended independently. The values
                listed here determine the indices where new
                segments start. The data will be split into
                len(bp)+1 segments. (Default: do not break input data)
            numpieces: Automatically determine breakpoints by splitting
                input data into roughly equal parts. This option, if provided,
                will override 'bp'. (Default: treat data as 1 piece).

        Output:
            detrended: a 1D array.
    """
    ymasked = np.ma.masked_array(ydata, mask=np.ma.getmaskarray(ydata))
    if xdata is None:
        xdata = np.ma.masked_array(np.arange(ydata.size), mask=np.ma.getmaskarray(ydata))
    detrended = ymasked.copy()
    
    if numpieces is None:
        edges = [0]+bp+[len(ydata)]
    else:
        # Determine indices to split at based on desired numbers of pieces
        isplit = np.linspace(0, len(ydata), numpieces+1, endpoint=1)
        edges = np.round(isplit).astype(int)
    for start, stop in zip(edges[:-1], edges[1:]):
        if not np.ma.count(ymasked[start:stop]):
            # No unmasked values, skip this segment.
            # It will be masked in the output anyway.
            continue
        x, poly_ydata = fit_poly(ymasked[start:stop], xdata[start:stop], order)
        detrended.data[start:stop] -= poly_ydata
    if np.ma.isMaskedArray(ydata):
        return detrended
    else:
        return detrended.data


def fit_poly(ydata, xdata, order=1):
    """Fit a polynomial to data using scipy.linalg.lstsq().
        
        Inputs:
            ydata: A 1D array to be detrended.
            xdata: A 1D array of x-values to use
            order: Order of polynomial to use (Default: 1)
        
        Outputs:
            x: An array of polynomial order+1 coefficients
            poly_ydata: A array of y-values of the polynomial evaluated 
                at the input xvalues.
    """
    # Convert inputs to masked arrays 
    # Note these arrays still reference the original data/arrays
    xmasked = np.ma.asarray(xdata)
    ymasked = np.ma.asarray(ydata)
    if not np.ma.count(ymasked):
        # No unmasked values!
        raise ValueError("Cannot fit polynomial to data. " \
                        "There are no unmasked values!")
    ycomp = ymasked.compressed()
    xcomp = xmasked.compressed()

    powers = np.arange(order+1)
 
    A = np.repeat(xcomp, order+1)
    A.shape = (xcomp.size, order+1)
    A = A**powers

    x, resids, rank, s = scipy.linalg.lstsq(A, ycomp)
    
    # Generate decompressed detrended array
    A = np.repeat(xmasked.data, order+1)
    A.shape = (len(xmasked.data), order+1)
    A = A**powers

    poly_ydata = np.dot(A, x).squeeze()
    
    return x, poly_ydata
