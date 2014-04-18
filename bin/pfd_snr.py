#!/usr/bin/env python
"""
Manually determine on-pulse region for the profile in a pfd file
    and calculate the SNR.

    Patrick Lazarus, May 14, 2013 (taken out of gridding.py).
"""

import sys
import argparse
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.integrate import trapz

import prepfold
import psr_utils

from pypulsar.utils.astro import protractor
from pypulsar.utils.astro import sextant
from pypulsar.utils import skytemp
debug = 1


class OnPulseError(Exception):
    pass

class Observation:
    """ Observation object
    """
    def __init__(self, pfd, sefd=None, verbose=True):
        """Return an observation object for the given pfd file.
    
        Inputs: 
            pfd: a pfd object
            sefd: the system-equivalent flux density of the observation (in Jy)
            verbose: if True, print extra information (Default: True)

        Output:
            obs: The Observation object            
        """
        self.sefd = sefd
        self.p = pfd 
        self.fn = self.p.pfd_filename
        print self.fn
        self.snr = None
        self.smean = None
        self.verbose = verbose
        self.notes = []
        
        self.p.dedisperse(doppler=True)
        self.p.adjust_period()
        prof = self.p.sumprof
        self.proflen = len(prof)
        self.nbin = len(prof)
        imax = np.argmax(prof)
        self.nrot = (imax-len(prof)/2) % len(prof)
        if self.verbose:
            print "Profile maximum at bin %d. Rotating by %d bins." % (imax, self.nrot)
        self.prof = np.asarray(psr_utils.rotate(prof, self.nrot))
        
        self.region_start = None
        self.region_start_line = None
        self.regions = []
        
        # Plot
        self.fig = plt.figure()
        self.ax = plt.gca()
        plt.plot(self.prof, 'k-', drawstyle='steps-post')
        
        # Set up triggers
        self.cid_mouseclick = self.fig.canvas.mpl_connect('button_press_event', \
                                                            self.mousepress)
        self.cid_pick = self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.cid_keypress = self.fig.canvas.mpl_connect('key_press_event', \
                                                            self.keypress)

    def mousepress(self, event):
        if event.inaxes and event.button==1:
            self.eventpress = event
            if self.region_start is None:
                # Starting a new on-pulse region
                xx =  int(event.xdata+0.5)
                self.region_start = xx
                self.region_start_line = plt.axvline(xx, c='k', ls='--')
            else:
                # Selected an on-pulse region
                if self.region_start > event.xdata:
                    xlo = 0
                    xhi = int(event.xdata+0.5)
                    self.regions.append(plt.axvspan(xlo, xhi, picker=True,
                                                    facecolor='b', alpha=0.5))
                    xlo = int(self.region_start+0.5)
                    xhi = self.nbin
                    self.regions.append(plt.axvspan(xlo, xhi, picker=True,
                                                    facecolor='b', alpha=0.5))
                else:
                    xlo = int(self.region_start+0.5)
                    xhi = int(event.xdata+0.5)
                    self.regions.append(plt.axvspan(xlo, xhi, picker=True,
                                                    facecolor='b', alpha=0.5))
                
                # Remove line
                self.region_start_line.remove()
                self.region_start = None
                self.region_start_line = None

            self.fig.canvas.draw()

    def onpick(self, event):
        if event.mouseevent.button==3:
            ind = self.regions.index(event.artist)
            self.regions.pop(ind)
            event.artist.remove()
            self.fig.canvas.draw()

    def calc_snr(self):
        ionpulse = np.zeros_like(self.prof, dtype=bool)
        for polygon in self.regions:
            xlo, xhi = polygon.get_xy()[0:3:2,0]
            ionpulse[xlo:xhi] = 1

        nbins_selected = ionpulse.sum()
        if nbins_selected < 1:
            warnings.warn("No on-pulse region selected!")
            return
        # Correct standard deviation for correlations between bins
        data_avg, data_var = self.p.stats.sum(axis=1).mean(axis=0)[1:3]
        nbin_eff = self.proflen*self.p.DOF_corr()
        std = np.sqrt(data_var*self.p.Nfolded/nbin_eff)
       
        # Calculate S/N using eq. 7.1 from Lorimer and Kramer
        offpulse = self.prof[~ionpulse]

        mean = offpulse.mean()
        scaled = self.prof-mean
        area = np.sum(scaled[ionpulse]) 
        profmax = np.max(scaled[ionpulse])
        weq = area/profmax
        self.snr = area/std/np.sqrt(weq)
        if self.verbose:
            print "Number of bins selected: %d (%f phase)" % \
                    (nbins_selected, nbins_selected/float(len(self.prof)))
            if debug:
                print "Equivalent width (bins):", weq
                print "Std-dev corrected for correlations between phase bins:", std
                print "Off-pulse mean:", mean
                print "Integral under the mean-subtracted on-pulse region:", \
                        area
            print "SNR:", self.snr

        if self.sefd is not None:
            npol = 2  # prepfold files only contain total-intensity
                      # (i.e. both polarisations summed)
            bw = self.p.chan_wid*self.p.numchan
            self.smean = self.snr*self.sefd/np.sqrt(npol*self.p.T*bw)*np.sqrt(weq/(len(self.prof)-weq))
            if self.verbose:
                print "Mean flux density (mJy):", self.smean
    
    def keypress(self, event):
        if event.key == ' ':
            self.calc_snr()
        elif event.key in ('q', 'Q'):
            self.calc_snr()
            plt.close(self.fig)
        elif event.key in ('r', 'R'):
            self.notes.append("RFI!")
        elif event.key in ('n', 'N'):
            self.notes.append("No detection!")


def main():
    for pfdfn in args.files:
        print pfdfn
        pfd = prepfold.pfd(pfdfn)
        if args.sefd is not None:
            sefd = args.sefd
        elif args.gain is not None and args.tsys is not None:
            fctr = 0.5*(pfd.hifreq+pfd.lofreq)
            glon, glat = sextant.equatorial_to_galactic(pfd.rastr, pfd.decstr, 
                                                        input='sexigesimal', 
                                                        output='deg')
            print "Galactic Coords: l=%g deg, b=%g deg" % (glon, glat)
            tsky = skytemp.get_skytemp(glon, glat, freq=fctr)[0]
            print "Sky temp at %g MHz: %g K" % (fctr, tsky)
            sefd = (args.tsys+tsky)/args.gain
        obs = Observation(pfd, sefd=sefd, verbose=True)
        plt.show()
        print " ".join(obs.notes)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate SNR from .pfd files.")
    parser.add_argument('files', nargs='*', \
                        help="Files to find SNR for.")
    parser.add_argument('--on-pulse', dest='on_pulse', nargs=2, type=float, \
                        help="On-pulse region. Two values should be provided, " \
                            "the starting phase and the ending phase.")
    parser.add_argument('--sefd', dest='sefd', type=float, \
                        help="The SEFD (in Jy) of the observing system. "
                             "(i.e. Tsys/Gain) Sky temperature will not be factored in.")
    parser.add_argument('--tsys', dest='tsys', type=float, \
                        help="The temperature of the observing system in K "
                             "(not including sky temperature).")
    parser.add_argument('--gain', dest='gain', type=float, \
                        help="The gain of the observing system in K/Jy.")
    args = parser.parse_args()

    if args.sefd is not None and (args.tsys is not None or args.gain is not None):
        raise ValueError("Gain and/or system temperature should not be " 
                         "provided if SEFD is given.")
    elif (args.tsys is not None and args.gain is None) or \
         (args.tsys is None and args.gain is not None):
        raise ValueError("Both gain and system temperature must be provided " 
                         "together.")
    main()
