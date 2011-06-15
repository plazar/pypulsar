#!/usr/bin/env python
"""autozap.py

Generate a zaplist by considering the median of multiple FFTs.

Patrick Lazarus, Nov. 17, 2010
"""

import glob
import sys
import os.path
import optparse

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import scipy.stats
import scipy.optimize

from pypulsar.utils import mydetrend
from pypulsar.formats import prestofft

# Define some constants
BLOCKSIZE = 10000
SMOOTHFACTOR = 10
MAXITER = 10

def get_ffts(fftfns):
    """Given a list of fftfilenames create PrestoFFT objects.
    """
    print "Number of .fft files found: %d" % len(fftfns)
   
    allpffts = [prestofft.PrestoFFT(fn, delayread=True, delayfreqs=True) \
		for fn in fftfns if not fn.endswith("7.fft")]
    print "Excluding %d FFTs of beam 7 data..." % (len(fftfns) - len(allpffts))
    
    bad_pffts = 0
    p1 = allpffts[0]
    p1size = os.path.getsize(p1.fftfn)
    pffts = []
    for pcurr in allpffts:
	if os.path.getsize(pcurr.fftfn) != p1size:
	    bad_pffts += 1
	else:
	    pffts.append(pcurr)

    if bad_pffts:
	print "Excluding %d FFTs of different size..." % bad_pffts
    print "Number of power spectra being considered: %d" % len(pffts)
    return pffts


def calc_percentile(pffts, percent=50):
    """Calculate percentile of power spectra from
	list of PrestoFFT objects.

        Inputs:
            pffts: a list of PrestoFFT objects
            percent: percent of FFTs with values smaller than value returned.
                    (A value between 0 and 100. Default = Calculate median)
    """
    p1 = pffts[0]
    pwrspec_size = p1.inf.N/2
    percentile = np.zeros(pwrspec_size)
    sys.stdout.write("Calculating percentile score (%.1f %%)... %5.1f %%" % (percent, 0))
    sys.stdout.flush()
    for block in np.arange(0,pwrspec_size, BLOCKSIZE):
	if block+BLOCKSIZE >= pwrspec_size:
	    blockend = pwrspec_size
	else:
	    blockend = block+BLOCKSIZE
	curr_block = []
        for pcurr in pffts:
            fft = pcurr.read_fft(count=(blockend-block))
            powers = np.abs(fft)**2
	    curr_block.append(powers)
	pwrspec_stack = np.array(curr_block)
	# Sort and reverse so values are descending
	percentile[block:blockend] = \
                scipy.stats.scoreatpercentile(pwrspec_stack, percent)
        sys.stdout.write("\rCalculating percentile score (%.1f %%)... %5.1f %%" % \
                            (percent, (100.0*blockend/pwrspec_size)))
        sys.stdout.flush()
    print "\rCalculating percentile score (%.1f %%)... Done       " % percent
    return percentile


def plot(freqs, spectrum, style="k-", **plotargs):
    """Plot the given spectrum.
    """
    # spectrum /= np.max(spectrum[(freqs>=0.125) & (freqs<1000)])

    subones = (freqs>=0.125) & (freqs<1) # 0.125 Hz is min freq searched in PALFA
    ones = (freqs>=1) & (freqs<10)
    tens = (freqs>=10) & (freqs<100)
    hundreds = (freqs>=100) & (freqs<1000)
    thousands = (freqs>=1000) & (freqs<10000)
    
    plt.rc('xtick', direction='out')
    plt.rc('xtick.major', size=4)
    plt.rc('xtick.minor', size=3)
    plt.rc(('xtick', 'ytick'), labelsize='x-small')
    
    fig = plt.gcf()
    if len(fig.get_axes())==5:
        print "Using existing axes"
        axsubones, axones, axtens, axhundreds, axthousands = fig.get_axes()
    else:
        fig.set_size_inches(10,8.5, forward=True)
        plt.subplots_adjust(hspace=0.3)
        print "Creating new axes"
        axsubones = plt.subplot(5,1,1) 
        axones = plt.subplot(5,1,2, sharey=axsubones)
        axtens = plt.subplot(5,1,3, sharey=axsubones)
        axhundreds = plt.subplot(5,1,4, sharey=axsubones)
        axthousands = plt.subplot(5,1,5, sharey=axsubones)

    # plot mean power spectrum (in five parts)
    plt.axes(axsubones)
    plt.plot(freqs[subones], spectrum[subones], style, **plotargs)
    plt.ylabel("Power")
    plt.xlim(0.1, 1)
    plt.xscale('log')
    axsubones.xaxis.set_ticks_position('bottom')

    plt.axes(axones)
    plt.plot(freqs[ones], spectrum[ones], style, **plotargs)
    plt.ylabel("Power")
    plt.xlim(1, 10)
    plt.xscale('log')
    axones.xaxis.set_ticks_position('bottom')
    
    plt.axes(axtens)
    plt.plot(freqs[tens], spectrum[tens], style, **plotargs)
    plt.ylabel("Power")
    plt.xlim(10, 100)
    plt.xscale('log')
    axtens.xaxis.set_ticks_position('bottom')
    
    plt.axes(axhundreds)
    plt.plot(freqs[hundreds], spectrum[hundreds], style, **plotargs)
    plt.ylabel("Power")
    plt.xlim(100, 1000)
    plt.xscale('log')
    axhundreds.xaxis.set_ticks_position('bottom')

    plt.axes(axthousands)
    plt.plot(freqs[thousands], spectrum[thousands], style, **plotargs)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power")
    plt.xlim(1000, 10000)
    plt.xscale('log')
    axhundreds.xaxis.set_ticks_position('bottom')


def gen_mask(freqs, powerspec, nsig=3.5):
    """Find peaks in the power spectrum. Return a boolean array
        of the same size as freqs where True values are peaks in 
        the power spectrum 'nsig'*sigma larger than the baseline
        of the spectrum.

        Inputs:
            freqs: an array of frequencies
            powerspec: a power spectrum (numpy array)
            nsig: Number of std deviations away from mean a bin has
                    to be in order to be zapped.

        Output:
            mask: a boolean array where True values are frequencies
                    to be zapped.
    """
    sys.stdout.write("Median filtering... ")
    sys.stdout.flush()
    filtered = scipy.signal.medfilt(powerspec, 101)
    print "Done"
    flattened = (powerspec-filtered)
    halfflat = flattened[flattened<0]
    halfflat.sort()
    cdfresids = lambda sigma: scipy.stats.norm(loc=0, scale=sigma).cdf(halfflat)-(np.arange(1, (halfflat.size+1))/(halfflat.size*2.0))
    sys.stdout.write("Optimizing... ")
    sys.stdout.flush()
    fitresult = scipy.optimize.leastsq(cdfresids, np.abs(np.array([halfflat[halfflat.size/2]])))
    sigma = fitresult[0]
    print "Done"
    mask = smooth(flattened, SMOOTHFACTOR) > (sigma*nsig)
    return mask
  

def hone_mask(freqs, powerspec, inmask, nsig):
    """Given a starting mask improve it.
        
        Inputs:
            freqs: an array of frequencies
            powerspec: a power spectrum (numpy array)
            inmask: A boolean array where True values are frequencies
                    to be zapped.
            nsig: Number of std deviations away from mean a bin has
                    to be in order to be zapped.

        Output:
            outmask: An improved mask.
    """
    outmask = np.zeros(powerspec.size, dtype='bool')
    # Break spectrum into pieces and detrend
    sys.stdout.write("Improving mask... %5.1f %%" % 0)
    sys.stdout.flush()
    for block in np.arange(0,powerspec.size, BLOCKSIZE):
	if block+BLOCKSIZE >= powerspec.size:
	    blockend = powerspec.size
	else:
	    blockend = block+BLOCKSIZE
        
        # Overlap blocks slightly so that smoothing doesn't cause
        # edges of blocks to be de-weighted.
        lo = 0
        hi = 0
        if block-SMOOTHFACTOR >= 0:
            lo = SMOOTHFACTOR
        if blockend+SMOOTHFACTOR < powerspec.size:
            hi = SMOOTHFACTOR
        
        spec_block = powerspec[block-lo:blockend+hi]
        freq_block = freqs[block-lo:blockend+hi]
        mask_block = inmask[block-lo:blockend+hi]
        detrended_block = mydetrend.detrend(np.log10(spec_block), \
                                        xdata=np.log10(freq_block), \
                                        mask=mask_block, order=2)
        unmasked_block = detrended_block[~mask_block]
        std_block = unmasked_block.std()
        smoothed_block = smooth(detrended_block,SMOOTHFACTOR)[lo:detrended_block.size-hi]
        outmask[block:blockend] = (smoothed_block > (std_block*nsig))
        sys.stdout.write("\rImproving mask... %5.1f %%" % (100.0*blockend/powerspec.size))
        sys.stdout.flush()
    print "\rImproving mask... Done       "
    return outmask


def smooth(data, smoothfactor=1, verbose=False):
    """Smooth data by convolving with tophat of width
	'smoothfactor' bins. The height of the tophat
	is chosen such that RMS = 1.

        If verbose is True, display a message.
    """
    if smoothfactor > 1:
        if verbose:
            print "Smoothing by %d bins..." % smoothfactor
	kernel = np.ones(smoothfactor, dtype='float32') / \
		    np.sqrt(smoothfactor)
	return scipy.signal.convolve(data, kernel, 'same', old_behavior=False)


def write_zaplist(zapfn, freqs, mask):
    """Write masked frequencies to a zaplist file.
        
        Inputs:
            zapfn: Filename of zaplist
            freqs: Array of frequency values
            mask: Array of boolean values where True represents
                    a frequency bin to zap.
    """
    zapfile = open(zapfn, 'w')
    zapfile.write("# This file was created automatically with autozap.py\n")
    zapfile.write("# Lines beginning with '#' are comments\n")
    zapfile.write("# Lines beginning with 'B' are barycentric freqs (i.e. PSR freqs)\n")
    zapfile.write("#                 Freq                 Width\n")
    zapfile.write("# --------------------  --------------------\n")
    badfreqs = np.ma.masked_array(freqs, mask=np.bitwise_not(mask))
    for s in np.ma.notmasked_contiguous(badfreqs):
        lofreq = freqs[s.start]
        hifreq = freqs[s.stop+1]
        width = (hifreq-lofreq)/2.0
        midfreq = (hifreq+lofreq)/2.0
        zapfile.write("  %20.15g  %20.15g\n" % (midfreq, width))
    zapfile.close()


def main():
    fftfns = glob.glob(options.glob) + sys.argv
    pffts = get_ffts(fftfns)

    # Calculate frequencies and power spectrum
    pffts[0].calcfreqs()
    freqs = pffts[0].freqs
    powerspec = calc_percentile(pffts, percent=options.percent)

    # Remove DC offset
    freqs = freqs[1:]
    powerspec = powerspec[1:]
    
    # Generate mask
    mask = gen_mask(freqs, powerspec, nsig=options.nsig)
    # Improve the mask
    for ii in range(MAXITER):
        newmask = hone_mask(freqs, powerspec, mask, options.nsig)
        if np.all(newmask==mask):
            print "Mask is stable."
            break
        else:
            mask = newmask

    # Write out the zaplist
    write_zaplist(options.outname+".zaplist", freqs, mask)
    
    # Plot
    plot(freqs, powerspec, style='r-', lw=0.25, zorder=-1)
    maskedspec = np.ma.masked_array(powerspec, mask=mask)
    plot(freqs, maskedspec, style='k-', lw=0.5, zorder=1)
    plt.suptitle("Percentile power spectrum (%.1f %%). " \
                 "Number of spectra combined: %d" % \
                    (options.percent, len(pffts)))
    plt.show()


if __name__ == '__main__':
    parser = optparse.OptionParser(usage="%prog [options] <fft fns>", \
                            version="v0.1 Patrick Lazarus, Nov. 18, 2010")
    parser.add_option('-g', '--glob', dest='glob', type='string', \
                        help="Glob expression referring to *.fft files. " \
                             "(The string must be properly quoted.)", \
                        default="")
    parser.add_option('--median', dest='percent', action='store_const', \
                        help="Compute median power spectrum. (This is " \
                             "equivalent to '{-p | --percent} 50').", \
                        const=50)
    parser.add_option('-p', '--percent', dest='percent', action='store', \
                        type='float', help="Computer percentile power " \
                             "spectrum. That is, the value at a given " \
                             "frequency is larger than PERCENT % of the " \
                             "input power spectra at the same frequency. " \
                             "(Default: 50, i.e. median power spectrum).", \
                         default=50)
    parser.add_option('-s', '--nsig', dest='nsig', type='float', \
                        help="Number of sigmas from mean to be considered " \
                             "an RFI spike in power spectrum. (Default: 3).", \
                        default=3)
    parser.add_option('-o', '--outname', dest='outname', \
                        help="Output filename's basename (i.e. no extension)", \
                        default="autozapped")
    (options, sys.argv) = parser.parse_args()
    main()
