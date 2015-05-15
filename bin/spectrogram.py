#!/usr/bin/env python

"""
Plot spectrogram (spin freq vs. time) for a datfile.

Patrick Lazarus, May 29, 2010
"""
import sys
import optparse

import matplotlib.pyplot as plt
import numpy as np

import datfile

def get_spectra(dat, time=1.0):
    """Break the timeseries in dat into blocks of length 'time' seconds.
        Compute the power spectrum of each block and return a 2D array
        of the spectra.

        Default is to produce a spectrum for each 1.0 second block.
    """
    samp_per_block = int(time/dat.infdata.dt)
    numspec = int(dat.infdata.N/samp_per_block)
    numcoeffs = samp_per_block/2+1

    spectra = np.empty((numspec, numcoeffs))
    samples = np.arange(numspec)*samp_per_block
    dat.rewind()
    for ii in np.arange(numspec):
        block = dat.read_Nsamples(samp_per_block)
        spectra[ii,:] = np.abs(np.fft.rfft(block))**2
    freqs = np.fft.helper.fftfreq(samp_per_block, dat.infdata.dt)
    freqs = freqs[freqs>=0]
    times = samples*dat.infdata.dt
    return spectra, times, freqs


def keypress(event):
    if event.key in ('q', 'Q'):
        sys.exit(0)


def main():
    datfn = sys.argv[0]
    dat = datfile.Datfile(datfn)
    spectra, times, freqs = get_spectra(dat, time=options.time)
    fig = plt.figure(figsize=(11,8.5))
    # Plot as image. Omit DC level (first coefficient in power spectra).
    spectrogram = spectra[:,1:]
    if options.log:
        spectrogram = np.log10(spectrogram)
    plt.imshow(spectrogram, aspect='auto', interpolation='bilinear', \
                extent=(freqs[1], freqs[-1], times[0], times[-1]))
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Time (s)")
    plt.title("Spectrogram of\n%s" % datfn)
    cb = plt.colorbar()
    if options.log:
        cb.set_label(r"log$_{10}$(Raw Power Spectrum Intensity)")
    else:
        cb.set_label("Raw Power Spectrum Intensity")
    plt.figtext(0.05, 0.025, "Integration time: %g s" % options.time, size='small')
    fig.canvas.mpl_connect('key_press_event', keypress)
    plt.show()


if __name__=='__main__':
    parser = optparse.OptionParser(prog="spectrogram.py", \
                        version="v0.9 Patrick Lazarus (May 29, 2010)")
    parser.add_option('-t', '--time', dest='time', type='float', \
                        help="Time duration (in seconds) of each block " \
                             "for which a power spectrum is calculated. " \
                             "(Default: 1 s.)", \
                        default=1.0)
    parser.add_option('-l', '--log', dest='log', action='store_true', \
                        help="Show logarithm colour scale for power " \
                             "spectrum intensity. (Default: Show linear " \
                             "colour scale.)", \
                        default=False)
    (options, sys.argv) = parser.parse_args()
    main()
