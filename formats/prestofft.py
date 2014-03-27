"""
prestofft.py

Read PRESTO .fft file.
"""

import sys
import warnings
import os.path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import infodata

"""
Define object.
frequencies list (generator?)
powers list
normalisation method? Keep copy of powers as well as raw powers?
inf file
plot method
"""

class PrestoFFT:
    def __init__(self, fftfn, inffn=None):
        """PrestoFFT object creator
            'fftfn' is filename of .fft file
            'inffn' is filename of .inf file
                If None the inffn will default to fftfn with .fft extension replaced by .inf
        """
        if not fftfn.endswith(".fft"):
            ValueError("FFT filename must end with '.fft'!")
        if not os.path.isfile(fftfn):
            ValueError("FFT file does not exist!\n\t(%s)" % fftfn)
        else:
            self.fftfn = fftfn
            self.fftfile = open(self.fftfn, "rb")

            if inffn is None:
                inffn = "%s.inf" % fftfn[:-4]
            else:
                inffn = inffn
            if not os.path.isfile(inffn):
                ValueError("Info file does not exist!\n\t(%s)" % inffn)

            self.inffn = inffn
            self.inf = infodata.infodata(inffn)

            freqs = np.fft.fftfreq(self.inf.N, self.inf.dt)
            self.freqs = freqs[freqs>=0]
            self.fft = self.read_fft()
            self.phases = np.angle(self.fft)

            self.normalisation = "raw"
            self.powers = np.abs(self.fft)**2

    def interpolate(self, r, m=32):
        """Interpolate the value of the FFT at real bin indices 'r'.
            Use 'm' nearest bins when interpolating.

            Inputs:
                r: real bin indices to interpolate the FFT at.
                m: The number of nearby bins to use. Must be an
                    even integer.

            Output:
                interpfft: FFT coefficients interpolated at 'r'.
        """
        if (m % 2) is not 0:
            raise ValueError("Input 'm' must be an even integer: %s" % str(m))
        round_r = np.round(r).astype('int')
        k = round_r[:,np.newaxis]+np.arange(-m/2, m/2+1)
        coefs = self.fft[k]
        expterm = np.exp(-1.0j*np.pi*(r[:,np.newaxis]-k))
        sincterm = np.sinc(np.pi*(r[:,np.newaxis]-k))
        interpfft = np.sum(coefs*expterm*sincterm, axis=1)
        return interpfft

    def read_fft(self, count=-1):
        """Read 'count' powers from .fft file and return them.
            power = real*real + imag*imag
        """
        fft = np.fromfile(self.fftfile, dtype=np.dtype('c8'), count=count)
        return fft

    def plot(self):
        """Plot the power spectrum.
        """
        plt.plot(self.freqs, self.powers, 'k-')
        plt.title(self.fftfn)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power")
        plt.show()

    def plot_3pane(self, zapfile=None):
        """Plot the power spectrum in 3 panes.
            If 'showzap' is True show PALFA zaplist.
        """
        ones = (self.freqs>=1) & (self.freqs<10)
        tens = (self.freqs>=10) & (self.freqs<100) 
        hundreds = (self.freqs>=100) & (self.freqs<1000)

        fig = plt.figure(figsize=(10,8))
        plt.subplots_adjust(hspace=0.25)
    
        # plot mean power spectrum (in three parts)
        axones = plt.subplot(3,1,1) 
        plt.plot(self.freqs[ones], self.powers[ones], 'k-', lw=0.5)
        plt.ylabel("Power") 
        plt.xscale('log')

        axtens = plt.subplot(3,1,2, sharey=axones)
        plt.plot(self.freqs[tens], self.powers[tens], 'k-', lw=0.5)
        plt.ylabel("Power")
        plt.xscale('log')

        axhundreds = plt.subplot(3,1,3, sharey=axones)
        plt.plot(self.freqs[hundreds], self.powers[hundreds], 'k-', lw=0.5)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power")
        plt.xscale('log')

        maxpwr = np.max(self.powers[(self.freqs>=1) & (self.freqs<1000)])
        axones.set_ylim(0, maxpwr*1.1)

        if zapfile is not None:
            # Plot regions that are zapped
            zaplist = np.loadtxt(zapfile)
            for freq, width in zaplist:
                for ax in (axones, axtens, axhundreds):
                    r = matplotlib.patches.Rectangle((freq-width/2.0, 0), width, maxpwr*1.1, \
                                                fill=True, fc='b', ec='none', \
                                                alpha=0.25, zorder=-1)
                    ax.add_patch(r)
            plt.figtext(0.025, 0.03, "Zaplist file: %s" % zapfile, size="xx-small")

        plt.suptitle("Power Spectrum (%s)" % self.fftfn)

        def close(event):
            if event.key in ('q','Q'):
                plt.close()
        fig.canvas.mpl_connect("key_press_event", close)

        plt.show()
 
