#!/usr/bin/env pythong

"""
Plot frequency vs. time (non-dedispersed) for a filterbank file.

This program is useful for determining if single pulses have
the appropriate dispersion relation.

Patrick Lazarus, April 27, 2010

TODO:
- Clean up code. Especially samples vs. downsampled samples.
- Add a summed profile. (DONE).
- Write function to do dedispersion. (DONE).
- Allow for combining channels by adding dedispersed adjacent channels.
- Allow user to move the trace (follows a dot on summed profile).
- Ensure order of frequencies in filterbank object is always high to low.
- Grab extra samples for summed profile (dedispersion).
"""
import os.path
import optparse
import sys

import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np
import scipy.signal

import fbobs
import psr_utils
import rfifind


def main():
    filfns = args # addtional argument is fileterbank files
    if options.debug:
        print "Input filterbank files:", filfns

    obs = fbobs.fbobs(filfns)
    # filfile.print_header()
    obslen = obs.obslen 

    # Determine start and end of interval (in seconds and samples)
    if options.start < 0:
        options.start = 0
    if (options.end is None) or (options.end > obslen):
        # set to end of filterbank file
        options.end = obslen

    reqstartsamp = int(options.start / obs.tsamp) # requested
    # Round down to a multiple of downsamp bins
    reqstartsamp = reqstartsamp - (reqstartsamp % options.downsamp)
    # Get extra bins for smoothing
    startsamp = reqstartsamp - options.width*options.downsamp

    reqendsamp = int(options.end / obs.tsamp) # requested
    # Round up to a multiple of downsamp bins
    reqendsamp = reqendsamp - (reqendsamp % options.downsamp) + options.downsamp
    # Get extra bins for smoothing
    endsamp = reqendsamp + options.width*options.downsamp
    if options.dm:
        # Get extra bins for dedispersion
        # Compute DM delays
        delay_seconds = psr_utils.delay_from_DM(options.dm, obs.frequencies)
        delay_seconds -= np.min(delay_seconds)
        delay_samples = delay_seconds/(options.downsamp*obs.tsamp)
        maxsamps = np.max(delay_samples*options.downsamp)
        maxsamps = int(np.round(float(maxsamps)/options.downsamp))*options.downsamp
        endsamp += maxsamps

    reqnumsaps = reqendsamp - reqstartsamp
    numsamps = endsamp - startsamp

    if options.debug:
        print "Requested start time: %s s (%d samples)" % \
                    (options.start, reqstartsamp)
        print "Actual start time: %s s (%d samples)" % \
                    (startsamp*obs.tsamp, startsamp)
        print "Requested end time: %s s (%d samples)" % \
                    (options.end, reqendsamp)
        print "Actual end time: %s s (%d samples)" % \
                    (endsamp*obs.tsamp, endsamp)

    # read data
    data = obs.get_sample_interval(startsamp, endsamp).astype('float32')
    obs.close_all()
    data.shape = (numsamps, obs.nchans)
   
    if options.mask is not None:
        if options.debug:
            print "Masking channels using %s" % options.mask
        # Mask channels
        mask = rfifind.rfifind(options.mask)
        maskchans = mask.mask_zap_chans
        maskchans = obs.nchans - 1 - np.array(list(maskchans))
        data = mask_channels(data, maskchans)

    # Modify data
    if options.downsamp > 1:
        if options.debug:
            print "Downsampling by %d bins" % options.downsamp
        data = downsample(data, factor=options.downsamp)
    if options.width > 1:
        if options.debug:
            print "Smoothing with boxcar %d bins wide" % options.width
        data = smooth(data, factor=options.width)[options.width:-options.width,:]
        startsamp += options.width*options.downsamp
        endsamp -= options.width*options.downsamp

    # plot data as an image
    fig = plt.figure()
    ax = plt.axes((0.15, 0.15, 0.8, 0.7))
    data_scaled = scale(data, indep=options.scaleindep)
    data_scaled = data_scaled[:-maxsamps/options.downsamp]
    endsamp -= maxsamps
    plt.imshow(data_scaled.transpose(), \
                    aspect='auto', cmap=matplotlib.cm.binary, interpolation='nearest', \
                    extent=(startsamp/options.downsamp, endsamp/options.downsamp, \
                                obs.frequencies[-1], obs.frequencies[0]))
    plt.xlabel("Sample")
    plt.ylabel("Observing frequency (MHz)")
    plt.suptitle("Frequency vs. Time")
    fig.text(0.05, 0.02, r"Start time: $\sim$ %s s, End time: $\sim$ %s s, " \
                "Downsampled: %d bins, Smoothed: %d bins" % \
                (options.start, options.end, options.downsamp, options.width), \
                ha="left", va="center", size="x-small")
    if options.dm is not None:
        xlim = plt.xlim()
        ylim = plt.ylim()
        
        # Plot dispersion delay trace
        plt.plot(startsamp/options.downsamp+delay_samples, obs.frequencies, \
                    'r-', lw=5, alpha=0.25)
        plt.xlim(xlim)
        plt.ylim(ylim)

        plt.axes((0.15, 0.85, 0.8, 0.1), sharex=ax)
        dedisp_prof = dedisperse(data, delay_samples)[:-maxsamps/options.downsamp]
        plt.plot(np.linspace(xlim[0],xlim[1],dedisp_prof.size), dedisp_prof, 'k-')
        plt.xticks([])
        plt.yticks([])
        plt.xlim(xlim)
    fig.canvas.mpl_connect('key_press_event', keypress)
    plt.show()


def keypress(event):
    if event.key in ('q', 'Q'):
        print "Closing..."
        plt.close()


def dedisperse(data, delays):
    """Given 'data', a 2D array containing frequency
        channels along axis 1 and time samples along
        axis 0, return a dedispersed timeseries at
        by delaying each channel by the appropriate
        number of samples, given in the list 'delays'. 
        
        NOTE: Channels are padded with zeros where
        necessary.
    """
    # Dedisperse interval and show profile
    dedisp_prof = np.zeros_like(data[:,0])
    for ii, delay in enumerate(delays):
        temp = data[delay:,ii]
        dedisp_prof[0:temp.size] += temp
    return dedisp_prof


def mask_channels(data, maskchans):
    """Given 'data', a 2D array containing frequency
       channels along axis 1 and time samples along
       axis 0, and a list of channels to mask,
       return the data with channels masked.
    """
    numsamps, numchans = data.shape
    for chan in maskchans:
        data[:,chan] = np.zeros(numsamps)
    return data


def downsample(data, factor=1):
    """Downsample each row in 'data', by 'factor' bins. 
        This is accomplished by summing adjacent bins. 

        NOTE: Excess bins are truncated.
    """
    if factor <= 1:
        return data
    oldlen = data.shape[0]
    newlen = int(oldlen/factor)
    numchans = data.shape[1]
    downsampled = np.zeros((newlen, numchans))

    for ii in np.arange(numchans):
        temp = data[:,ii][:newlen*factor].copy()
        temp.shape = (newlen, factor)
        downsampled[:,ii] = temp.sum(axis=1)
    return downsampled


def smooth(data, factor=1):
    """Smooth each row in 'data', by 'factor' bins. 
        This is accomplished by convolving with a boxcar
        of width 'factor' bins. 

        NOTE: Data is padded with zeros, so bins at
              beginning and end of interval are supressed.
    """
    if factor <= 1:
        return data
    numchans = data.shape[1]
    kernel =  np.ones(factor, dtype='float32') / np.sqrt(factor)
    
    for ii in np.arange(numchans):
        data[:,ii] = scipy.signal.convolve(data[:,ii], kernel, 'same')
    return data

def scale(data, indep=False):
    """ Given 'data', a 2D array containing frequency 
        channels along axis 1 and time samples along 
        axis 0, scale 'data' and return the result.

        If 'indep' is True, each channel will be scaled
        independently.
    """
    numchans = data.shape[1]
    for ii in np.arange(numchans):
        data[:,ii] = (data[:,ii] - np.min(data[:,ii]))
        if indep:
            max = np.max(data[:,ii])
            if max != 0:
                data[:,ii] /= max
    if not indep:
        data /= np.max(data)
        
    return data


if __name__ == '__main__':
    parser = optparse.OptionParser(usage="%prog [options] FILTERBANK_FILE", \
                    description="Given a SIGPROG format filterbank file " \
                        "plot frequency vs. time (non-dedispersed) for " \
                        "a particular interval to verify if single pulses " \
                        "have the appropriate dispersion delay.",
                    version="%prof v0.9", prog="freq_time.py")
    parser.add_option('--debug', dest='debug', action='store_true', \
                    help="Display debugging information. (Default: don't " \
                         "display debugging information).", \
                    default=False)
    parser.add_option('--downsamp', dest='downsamp', type='int', \
                    help="Factor to downsample data by. (Default: 1).", \
                    default=1)
    parser.add_option('-w', '--width', dest='width', type='int', \
                    help="Width, in samples, of boxcar to convolve with. " \
                         "(Default: 1).", \
                    default=1)
    parser.add_option('--dm', dest='dm', type='float', \
                    help="DM, in cm^-3 pc, used for dispersion delay trace. " \
                         "(Default: No trace).", \
                    default=None)
    parser.add_option('-s', '--start', dest='start', type='float', \
                    help="Start of interval to plot, in seconds elapsed " \
                         "since start of observation. (Default: 0 s).", \
                    default=0.0)
    parser.add_option('-e', '--end', dest='end', type='float', \
                    help="End of interval to plot, in second elapsed "\
                         "since start of observation. (Default: end of file).", \
                    default=None)
    parser.add_option('--mask', dest='mask', type='string', \
                    help="Mask file used to mask channels. (Default: No Mask).", \
                    default=None)
    parser.add_option('--scaleindep', dest='scaleindep', action='store_true', \
                    help="If this flag is set scale each channel independently. " \
                         "(Default: Scale using global maximum.)", \
                    default=False)
    options, args = parser.parse_args()
    main()
