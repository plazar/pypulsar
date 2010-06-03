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
    def __init__(self, fftfn, inffn=None, delayfreqs=True, delayread=False):
        """PrestoFFT object creator
            'fftfn' is filename of .fft file
            'inffn' is filename of .inf file
                If None the inffn will default to fftfn with .fft extension replaced by .inf
            'delayfreqs', if True the frequencies will not be calculated upon object creation.
            'delayread', if True the fftfile will not be read upon object creation
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

            self.freqs = None
            self.rawpowers = None
            if not delayfreqs:
                # Calculate frequencies immediately
                self.calcfreqs()
            if not delayread:
                # Read all data immediately
                self.rawpowers = self.read_fft()
            
            self.normalisation = "raw"
            self.powers = self.rawpowers

    def calcfreqs(self):
        """Calculate frequencies corresponding to power
            spectrum bins, and store as attribute of self.
        """
        self.freqs = np.arange(0, self.inf.N/2)/(self.inf.N*self.inf.dt)

    def read_fft(self, count=-1):
        """Read 'count' powers from .fft file and return them.
            power = real*real + imag*imag
        """
        fft = np.fromfile(self.fftfile, dtype=np.dtype('c8'), count=count)
        powers = np.abs(fft)**2
        return powers

    def plot(self):
        """Plot the power spectrum.
        """
        plt.plot(self.freqs, self.powers, 'k-')
        plt.title(self.fftfn)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Normalised Power")
        plt.show()

    def plot_3pane(self, showzap=False):
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

        if showzap:
            # Plot regions that are zapped
            zaplist = np.loadtxt("/homes/borgii/alfa/svn/workingcopy_PL/PALFA/miscellaneous/PALFA.zaplist")
            for freq, width in zaplist:
                for ax in (axones, axtens, axhundreds):
                    r = matplotlib.patches.Rectangle((freq-width/2.0, 0), width, maxpwr*1.1, \
                                                fill=True, fc='b', ec='none', \
                                                alpha=0.25, zorder=-1)
                    ax.add_patch(r)

        plt.suptitle("Power Spectrum (%s)" % self.fftfn)

        def close(event):
            if event.key in ('q','Q'):
                plt.close()
        fig.canvas.mpl_connect("key_press_event", close)

        plt.show()
 
