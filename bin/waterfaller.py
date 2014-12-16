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

SWEEP_STYLES = ['r-', 'b-', 'g-', 'm-', 'c-']

def get_filterbank_data(fil, startsamp, N):
    """Return 2D array of data from filterbank file.

        Inputs:
            fil: Filterbank object
            startsamp: Starting sample
            N: number of samples to read

        Output:
            data: 2D numpy array
    """
    fil.seek_to_sample(startsamp)
    data = fil.read_Nsamples(N)
    data.shape = (N, fil.header['nchans'])
    return spectra.Spectra(fil.frequencies, fil.tsamp, data.T, \
                            starttime=fil.tsamp*startsamp, dm=0)


def get_psrfits_data(pfits, startsamp, N):
    """Return 2D array of data from PSRFITS file.

        Inputs:
            pfits: PsrfitsFile object.
            startsamp, Starting sample
            N: number of samples to read

        Output:
            data: 2D numpy array
    """
    # Calculate starting subint and ending subint
    startsub = int(startsamp/pfits.nsamp_per_subint)
    ntrim_start = startsamp - (startsub*pfits.nsamp_per_subint)
    endsub = int((startsamp+N)/pfits.nsamp_per_subint)
    ntrim_end = ((endsub+1)*pfits.nsamp_per_subint) - (startsamp+N)
    # Read data
    data = pfits.read_subint(startsub)
    for isub in xrange(startsub+1, endsub+1):
        tmp = pfits.read_subint(isub)
        data = np.concatenate([data, tmp])
    # Truncate data to desired interval
    data = data[ntrim_start:]
    # Trim from end of array only if necessary
    if ntrim_end > 0:
        data = data[:-ntrim_end]
    return spectra.Spectra(pfits.freqs, pfits.tsamp, data, \
                           starttime=pfits.tsamp*startsamp, dm=0)


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
    blocknums = np.floor(sampnums/rfimask.ptsperint).astype('int')
    mask = np.zeros((N, rfimask.nchan), dtype='bool')
    for blocknum in np.unique(blocknums):
        blockmask = np.zeros_like(mask[blocknums==blocknum])
        blockmask[:,rfimask.mask_zap_chans_per_int[blocknum]] = True
        mask[blocknums==blocknum] = blockmask
    return mask.T
        

def main():
    fn = args[0]
    if fn.endswith(".fil"):
        # Filterbank file
        filetype = "filterbank"
        rawdatafile = filterbank.filterbank(fn)
        get_data = get_filterbank_data
    elif fn.endswith(".fits"):
        # PSRFITS file
        filetype = "psrfits"
        rawdatafile = psrfits.PsrfitsFile(fn)
        get_data = get_psrfits_data

    # Read data
    start_bin = np.round(options.start/fil.tsamp).astype('int')
    if options.nbins is None:
        nbins = np.round(options.duration/fil.tsamp).astype('int')
    else:
        nbins = options.nbins
    if options.maskfile is not None:
        rfimask = rfifind.rfifind(options.maskfile) 
        mask = get_mask(rfimask, start_bin, nbins)
        data = get_data(rawdatafile, start_bin, nbins)
        # Mask data
        data = data.masked(mask, maskval='median-mid80')
    else:
        data = get_data(rawdatafile, start_bin, nbins)

    # Subband data
    if (options.nsub is not None) and (options.subdm is not None):
        data.subband(options.nsub, options.subdm, padval='mean')

    # Dedisperse
    if options.dm:
        data.dedisperse(options.dm, padval='mean')

    # Downsample
    data.downsample(options.downsamp)

    # scale data
    data = data.scaled(options.scaleindep)
    
    # Smooth
    if options.width_bins > 1:
        data.smooth(options.width_bins, padval='mean')
    
    # Ploting it up
    fig = plt.figure()
    fig.canvas.set_window_title("Frequency vs. Time")
    ax = plt.axes((0.15, 0.15, 0.8, 0.7))
    plt.imshow(data.data, aspect='auto', \
                cmap=matplotlib.cm.cmap_d[options.cmap], \
                interpolation='nearest', origin='upper', \
                extent=(data.starttime, data.starttime+data.numspectra*data.dt, \
                        data.freqs.min(), data.freqs.max()))
    if options.show_cb:
        cb = plt.colorbar()
        cb.set_label("Scaled signal intensity (arbitrary units)")

    plt.axis('tight')

    # Sweeping it up
    for ii, sweep_dm in enumerate(options.sweep_dms):
        ddm = sweep_dm-data.dm
        delays = psr_utils.delay_from_DM(ddm, data.freqs)
        delays -= delays.min()
        
        if options.sweep_posns is None:
            sweep_posn = 0.0
        elif len(options.sweep_posns) == 1:
            sweep_posn = options.sweep_posns[0]
        else:
            sweep_posn = options.sweep_posns[ii]
        sweepstart = data.dt*data.numspectra*sweep_posn+data.starttime
        sty = SWEEP_STYLES[ii%len(SWEEP_STYLES)]
        plt.plot(delays+sweepstart, data.freqs, sty, lw=4, alpha=0.5)

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
    parser.add_option('--sweep-dm', dest='sweep_dms', type='float', \
                        action='append', \
                        help="Show the frequency sweep using this DM. " \
                                "(Default: Don't show sweep)", default=[])
    parser.add_option('--sweep-posn', dest='sweep_posns', type='float', \
                        action='append', \
                        help="Show the frequency sweep at this position. " \
                                "The position refers to the high-frequency " \
                                "edge of the plot. Also, the position should " \
                                "be a number between 0 and 1, where 0 is the " \
                                "left edge of the plot. "
                                "(Default: 0)", default=None)
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
