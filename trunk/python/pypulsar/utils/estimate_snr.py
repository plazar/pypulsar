"""
estimate_snr.py

A module for estimating signal-to-noise ratios.

Patrick Lazarus, November 1st, 2009
"""

import sys
import types
import numpy as np
from scipy import special

import skytemp

verbose = True

class SnrEstimator:
    """
    Signal-to-noise estimator object.
    """
    def __init__(self, freq, bw, numpol, gain, systemp, fwhm):
        """SnrEstimator object creator.

        freq: Centre frequency of observation (in MHz).
        bw: Bandwidth of observation (in MHz).
        numpol: Number of polarizations summed.
        gain: Telescope's gain (in Jy/K).
               * Can be a scalar, or a function of zenth and azimuth.
        systemp: Receiver's temperature (in K).
               * Can be a scalar, or a function of zenth and azimuth.
        fwhm: Full-width half-max of telescope's beam (in arcmin).
               * Can be a scalar, or a function of zenth and azimuth.
        """
        # Initialise values
        self.freq = freq
        self.bw = bw
        self.numpol = numpol
        
        if type(gain) is types.FunctionType:
            self.gain = gain
        else:
            self.gain = lambda za=0, az=0: gain
        
        if type(systemp) is types.FunctionType:
            self.systemp = systemp
        else:
            self.systemp = lambda za=0, az=0: systemp
        
        if type(fwhm) is types.FunctionType:
            self.fwhm = fwhm
        else:
            self.fwhm = lambda za=0, az=0: fwhm
      
        self.beam_profile = airy_pattern
      
    def estimate_snr(self, za, az, Smean, Sfreq, time, angsep, period, \
                        w50=None, Serror=None, l=None, b=None):
        """Estimate signal-to-noise of a known pulsar.
        Returns signal-to-noise and error on signal-to-noise.

        za: zenith angle (in deg).
        az: azimuth (in deg).
        Smean: mean flux density (in mJy).
        Sfreq: frequency at which Smean is given (in MHz).
        time: length of observation (in s).
        angsep: offset between beam centre and pulsar (in arcmin).
        period: pulsar's period (in s).
        w50: pulse width at 50% max intensity (in s).
                * default is 5% of pulse period, if not provided
        Serror: error on flux density (in mJy).
        l: galactic longitude of pulsar (in deg).
        b: galactic latitude of pulsar (in deg).
        """
        if w50 is None:
            w50 = 0.05*period

        if self.freq != Sfreq:
            Smean, Serror = change_freq(Smean, error=Serror, \
                                            oldfreq=Sfreq, newfreq=self.freq)
            Sfreq = self.freq
        
        if l is not None and b is not None:
            Tsky = skytemp.get_skytemp(l, b, freq=self.freq)
        else:
            if verbose:
                sys.stderr.write("No coords provided. Setting Tsky to 0 K.")
            Tsky = 0
        temp = self.systemp(za, az)+Tsky

        k = self.gain(za, az)*airy_pattern(self.fwhm(za, az), angsep) * \
            np.sqrt(self.numpol*time*self.bw)/temp * \
            np.sqrt((period-w50)/w50)

        Smean = np.atleast_1d(Smean)
        Serror = np.atleast_1d(Serror)
        
        snr = Smean*k
        snrerror = Serror*k
        snrerror[Serror is None or Serror==0] = np.nan
        return snr, snrerror


def airy_pattern(fwhm, x):
    """Return value of Airy pattern at x.
        Airy pattern is normalised so Airy(0)=1

        fwhm: full-width half-max of Airy pattern.
        x: distance from origin to evaluate
    """
    x = np.atleast_1d(x)
    # half-max of Airy pattern's central peak occurs at 1.61633
    scaled_x = x/fwhm*(2.0*1.61633)
    airy = np.atleast_1d((2*special.j1(scaled_x)/scaled_x)**2)
    airy[x==0] = 1
    return airy


def change_freq(S, error, oldfreq, newfreq, index=-1.8):
    """Change frequency.

    Return flux density and error for a new frequency.

    S: flux desnity (in mJy).
    error: error on flux density (in mJy).
    oldfreq: frequency of input S (in MHz).
    newfreq: frequency of output S (in MHz).
    [index]: spectral index to use.
    """
    k = (newfreq/oldfreq)**index
    newS = S*k
    if error is not None:
        newerror = error*k
    else:
        newerror = None
    return newS, newerror
