import numpy as np
import matplotlib.pyplot as plt

def hist(xx, bins, tot=None, bottom=None, 
         *args, **kwargs):
    """
        Inputs:
            xx: The data values to bin and plot a histogram for.
            bins: The number of bins to use, or the bin edges. 
                (This is passed directly to np.histogram)
            tot: The total number of points to use for scaling. 
                (Default: scale by total number of values in xx)
            bottom: The bottom edge of the histogram. 
                (Useful for stacked histograms)
    
            ***Additional arguments are passed directly to plt.fill

        Outputs:
            counts: The number of elements per bin.
            edges: The bin edges.
    """
    if tot is None:
        tot = float(len(xx))
    else:
        tot = float(tot)
    counts, edges = np.histogram(xx, bins=bins)
    counts = counts/tot
    if bottom is not None:
        counts += bottom
    x = bins.repeat(2)
    y = np.zeros(len(bins)*2)
    y[1:-1] = counts.repeat(2)
    plt.fill(x, y, *args, **kwargs)
    return counts, edges

