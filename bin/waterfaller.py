#!/usr/bin/env python

"""
waterfaller.py

Make waterfall plots to show frequency sweep of a single pulse.
Reads SIGPROC filterbank format.

Patrick Lazarus - Aug. 19, 2011
"""

import sys
import optparse

import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np

import psr_utils
import rfifind

from pypulsar.formats import filterbank
from pypulsar.formats import spectra

def get_data(fil, startsamp, N, mask=np.ma.nomask):
    """Return 2D array of data from filterbank file.

        Inputs:
            fil: Filterbank object
            startsamp: Starting sample
            N: number of samples to read
            mask: mask to use to remove bad values from data
                (Default: No mask)

        Output:
            data: 2D numpy array
    """
    fil.seek_to_sample(startsamp)
    data = fil.read_Nsamples(N)
    data.shape = (N, fil.header['nchans'])
    return spectra.Spectra(fil.frequencies, fil.tsamp, data.T, \
                            starttime=options.start, dm=0, mask=mask)


def get_mask(rfimask, startsamp, N):
    """Return an array of boolean values to act as a mask
        for a Spectra object.

        Inputs:
            rfimask: An rfifind.rfinfind object
            startsamp: Starting sample
            N: number of samples to read

        Output:
            mask: 2D numpy array of boolean values. 
                True represents an element that should be masked.
    """
    sampnums = np.arange(startsamp, startsamp+N)
    intnums = np.floor(sampnums/rfimask.ptsperint).astype('int')
    mask = np.zeros((N, rfimask.nchan), dtype='bool')
    for intnum in np.unique(intnums):
        mask[intnums==intnum][:,rfimask.mask_zap_chans_per_int[intnum]] = True
    return mask
        

def main():
    # Read filterbank file
    fil = filterbank.filterbank(args[0])
    
    # Read data
    start_bin = np.round(options.start/fil.tsamp).astype('int')
    if options.nbins is None:
        nbins = np.round(options.duration/fil.tsamp).astype('int')
    else:
        nbins = options.nbins
    if options.maskfile is not None:
        rfimask = rfifind.rfifind(options.maskfile) 
        mask = get_mask(rfimask, start_bin, nbins)
        data = get_data(fil, start_bin, nbins, mask=mask)
    else:
        data = get_data(fil, start_bin, nbins)

    # Subband data
    if (options.nsub is not None) and (options.subdm is not None):
        data.subband(options.nsub, options.subdm, padval='mean')
    print type(data.data)

    # Dedisperse
    if options.dm:
        data.dedisperse(options.dm, padval='mean')
    print type(data.data)

    # Downsample
    data.downsample(options.downsamp)
    print type(data.data)

    # scale data
    data = data.scaled(options.scaleindep)
    print type(data.data)
    
    # Smooth
    if options.width_bins > 1:
        data.smooth(options.width_bins, padval='mean')
    print type(data.data)
    
    # Ploting it up
    fig = plt.figure()
    fig.canvas.set_window_title("Frequency vs. Time")
    ax = plt.axes((0.15, 0.15, 0.8, 0.7))
    plt.imshow(data.data.mask, aspect='auto', \
                cmap=matplotlib.cm.cmap_d[options.cmap], \
                interpolation='nearest', origin='upper', \
                extent=(data.starttime, data.starttime+data.numspectra*data.dt, \
                        data.freqs.min(), data.freqs.max()))
    if options.show_cb:
        cb = plt.colorbar()
        cb.set_label("Scaled signal intensity (arbitrary units)")

    plt.axis('tight')

    # Sweeping it up
    if options.sweep_dm is not None:
        ddm = options.sweep_dm-data.dm
        delays = psr_utils.delay_from_DM(ddm, data.freqs)
        delays -= delays.min()
        
        sweepstart = data.dt*data.numspectra*options.sweep_posn+data.starttime
        plt.plot(delays+sweepstart, data.freqs, 'r-', lw=4, alpha=0.5)

    # Dressing it up
    plt.xlabel("Time")
    plt.ylabel("Observing frequency (MHz)")
    plt.suptitle("Frequency vs. Time")
    fig.canvas.mpl_connect('key_press_event', \
            lambda ev: (ev.key in ('q','Q') and plt.close(fig)))
    plt.show()


if __name__=='__main__':
    parser = optparse.OptionParser(prog="waterfaller.py", \
                        version="v0.9 Patrick Lazarus (Aug. 19, 2011)", \
                        usage="%prog [OPTIONS] INFILE", \
                        description="Create a waterfall plot to show the " \
                                    "frequency sweep of a single pulse " \
                                    "in SIGPROC filterbank data.")
    parser.add_option('--subdm', dest='subdm', type='float', \
                        help="DM to use when subbanding. (Default: " \
                                "same as --dm)", default=None)
    parser.add_option('-s', '--nsub', dest='nsub', type='int', \
                        help="Number of subbands to use. Must be a factor " \
                                "of number of channels. (Default: " \
                                "number of channels)", default=None)
    parser.add_option('-d', '--dm', dest='dm', type='float', \
                        help="DM to use when dedispersing data for plot. " \
                                "(Default: 0 pc/cm^3)", default=0.0)
    parser.add_option('-T', '--start-time', dest='start', type='float', \
                        help="Time into observation (in seconds) at which " \
                                "to start plot.")
    parser.add_option('-t', '--duration', dest='duration', type='float', \
                        help="Duration (in seconds) of plot.")
    parser.add_option('-n', '--nbins', dest='nbins', type='int', \
                        help="Number of time bins to plot. This option takes " \
                                "precedence over -t/--duration if both are " \
                                "provided.")
    #parser.add_option('-w', '--width-time', dest='width_time', type='float', \
    #                    help="Smooth each channel/subband with boxcar " \
    #                            "of this duration (in seconds). " \
    #                            "(Default: Don't smooth)")
    parser.add_option('--width-bins', dest='width_bins', type='int', \
                        help="Smooth each channel/subband with a boxcar " \
                                "this many bins wide. (Default: Don't smooth)", \
                        default=1)
    parser.add_option('--sweep-dm', dest='sweep_dm', type='float', \
                        help="Show the frequency sweep using this DM. " \
                                "(Default: Don't show sweep)", default=None)
    parser.add_option('--sweep-posn', dest='sweep_posn', type='float', \
                        help="Show the frequency sweep at this position. " \
                                "The position refers to the high-frequency " \
                                "edge of the plot. Also, the position should " \
                                "be a number between 0 and 1, where 0 is the " \
                                "left edge of the plot. "
                                "(Default: 0)", default=0.0)
    parser.add_option('--downsamp', dest='downsamp', type='int', \
                        help="Factor to downsample data by. (Default: 1).", \
                        default=1)
    parser.add_option('--mask', dest='maskfile', type='string', \
                        help="Mask file produced by rfifind. (Default: No Mask).", \
                        default=None)
    parser.add_option('--scaleindep', dest='scaleindep', action='store_true', \
                        help="If this flag is set scale each channel " \
                                "independently. (Default: Scale using " \
                                "global maximum.)", \
                        default=False)
    parser.add_option('--show-colour-bar', dest='show_cb', action='store_true', \
                        help="If this flag is set show a colour bar. " \
                                "(Default: No colour bar.)", \
                        default=False)
    parser.add_option('--colour-map', dest='cmap', \
                        help="The name of a valid matplotlib colour map." \
                                "(Default: gist_yarg.)", \
                        default='gist_yarg')
    options, args = parser.parse_args()
    
    if not hasattr(options, 'start'):
        raise ValueError("Start time (-T/--start-time) " \
                            "must be given on command line!")
    if (not hasattr(options, 'duration')) and (not hasattr(options, 'nbins')):
        raise ValueError("One of duration (-t/--duration) " \
                            "and num bins (-n/--nbins)" \
                            "must be given on command line!")
    if options.subdm is None:
        options.subdm = options.dm
    main()
