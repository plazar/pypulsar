#!/usr/bin/env python
"""
Manually determine on-pulse region for the profile in a pfd file
    and calculate the SNR.

    Patrick Lazarus, May 14, 2013 (taken out of gridding.py).
"""

import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.integrate import trapz

import prepfold
import psr_utils

from pypulsar.utils.astro import protractor
debug = 1


class OnPulseError(Exception):
    pass

class Observation:
    """ Observation object
    """
    def __init__(self, pfdfn, sefd=None):
        """Return an observation object for the given pfd file.
    
        Inputs: 
            pfdfn: a pfd filename
            sefd: the system-equivalent flux density of the observation (in Jy)

        Output:
            obs: The Observation object            
        """
        self.fn = pfdfn
        self.sefd = sefd
        self.p = prepfold.pfd(pfdfn)
        # RA in arcmin
        self.ra = protractor.convert(self.p.rastr, 'hmsstr', 'deg')*60
        # Dec in arcmin
        self.dec = protractor.convert(self.p.decstr, 'dmsstr', 'deg')*60
        self.snr = None
        self.width = None
        
        prof = self.p.bestprof.profile
        imax = np.argmax(prof)
        self.nrot = (imax-len(prof)/2) % len(prof)
        print "Profile maximum at bin %d. Rotating by %d bins." % (imax, self.nrot)
        self.prof = psr_utils.rotate(prof, self.nrot)

    def get_snr(self):
        """
        Return signal-to-noise ratio for the pfd object.

        Input: pfd - a pfd object.
        """
        self.fig = plt.figure()
        self.ax = plt.gca()
        plt.plot(self.prof, 'k-', drawstyle='steps-post')
        self.cid_keypress = self.fig.canvas.mpl_connect('key_press_event', \
                                                            self.keypress)
        self.cid_mouseclick = self.fig.canvas.mpl_connect('button_press_event', \
                                                            self.mousepress)
        self.rectprops = dict(facecolor='blue', edgecolor='blue', linewidth=2, \
                                alpha=0.5, fill=True)
        ymin, ymax = self.ax.get_ylim()
        self.to_draw = Rectangle((0,ymin-1e10), 0, ymax-ymin+2e10, visible=False, \
                                    **self.rectprops)
        self.ax.add_patch(self.to_draw)
        self.eventpress = None

        plt.show()

    def mousepress(self, event):
        if event.inaxes and event.button==1:
            self.mousepressed = True
            self.cid_mousemove = self.fig.canvas.mpl_connect('motion_notify_event', \
                                                                self.mousemove)
            self.cid_mouserelease = self.fig.canvas.mpl_connect('button_release_event', \
                                                                self.mouserelease)
            self.eventpress = event
            self.to_draw.set_visible(True)
            self.fig.canvas.draw()

    def mousemove(self, event):
        # update rectangle, draw
        if event.inaxes:
            xmin = self.eventpress.xdata
            xmax = event.xdata
            # Switch xmin/xmax
            if xmin > xmax: 
                xmin, xmax = xmax, xmin
            if xmin < 0:
                xmin = 0
            if xmax > len(self.prof)-1:
                xmax = len(self.prof)-1
            xmin = int(np.round(xmin))
            xmax = int(np.round(xmax))
            self.to_draw.set_x(xmin)
            self.to_draw.set_width(xmax-xmin)
            self.fig.canvas.draw_idle()
        
    def mouserelease(self, event):
        # remove rectangle (flash in red first for effect)
        self.to_draw.set_facecolor('red')
        self.to_draw.set_edgecolor('red')
        self.fig.canvas.draw()
        self.to_draw.set_visible(False)
        self.fig.canvas.draw()
        self.to_draw.set_facecolor('blue')
        self.to_draw.set_edgecolor('blue')
            
        # compute snr and pulse width
        xmin = self.eventpress.xdata
        xmax = event.xdata
        # Switch xmin/xmax
        if xmin > xmax: 
            xmin, xmax = xmax, xmin
        if xmin < 0:
            xmin = 0
        if xmax > len(self.prof)-1:
            xmax = len(self.prof)-1
        xmin = int(np.round(xmin))
        xmax = int(np.round(xmax))
        self.onpulse_range = (xmin, xmax)
        print self.onpulse_range
        self.calc_snr()
        self.disconnect_triggers()

    def calc_snr(self):
        if self.onpulse_range is None:
            raise OnPulseError("No on-pulse region selected!")
        # Correct standard deviation for correlations between bins
        nbin_eff = self.p.bestprof.proflen*self.p.DOF_corr()
        std = self.p.bestprof.data_std*np.sqrt(self.p.bestprof.N/nbin_eff)
       
        xmin, xmax = self.onpulse_range
        # Calculate S/N using eq. 7.1 from Lorimer and Kramer
        offpulse = np.concatenate((self.prof[:xmin], \
                                    self.prof[xmax:]))
        mean = offpulse.mean()
        scaled = self.prof-mean
        #print xmin, xmax
        #print scaled[xmin:xmax]
        #plt.figure()
        #plt.plot(scaled, drawstyle='steps-post')
        #plt.show()
        area = np.sum(scaled[xmin:xmax]) 
        profmax = np.max(scaled[xmin:xmax])
        weq = area/profmax
        self.width = xmax-xmin
        self.snr = area/std/np.sqrt(weq)
        print "Width selected: %d bins, %f phase" % \
                (self.width, self.width/float(len(self.prof)))
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
            smean = self.snr*self.sefd/np.sqrt(npol*self.p.T*bw)*np.sqrt(weq/(len(self.prof)-weq))
            print "Mean flux density (mJy):", smean

    def disconnect_triggers(self):
        # disconnect mousemove, mouserelease and reset state
        self.eventpress = None
        self.mousepressed = False
        self.fig.canvas.mpl_disconnect(self.cid_mousemove)
        self.fig.canvas.mpl_disconnect(self.cid_mouserelease)

    def keypress(self, event):
        if event.key == 'enter' and self.snr is not None and \
                self.width is not None:
            plt.close(event.canvas.figure.number)


def main():
    for pfdfn in args.files:
        print pfdfn
        obs = Observation(pfdfn, sefd=args.sefd)
        if args.on_pulse is None:
            obs.get_snr()
        else:
            print args.on_pulse
            start_bin = int(args.on_pulse[0]*len(obs.prof)-obs.nrot + 0.5) % len(obs.prof)
            end_bin = int(args.on_pulse[1]*len(obs.prof)-obs.nrot+0.5) % len(obs.prof)
            obs.onpulse_range = (start_bin, end_bin)
            print "Rotated on-pulse range (bins): %d - %d" % (start_bin, end_bin)
            obs.calc_snr()
        print "%s \tSNR: %s" % (obs.fn, obs.snr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate SNR from .pfd files.")
    parser.add_argument('files', nargs='*', \
                        help="Files to find SNR for.")
    parser.add_argument('--on-pulse', dest='on_pulse', nargs=2, type=float, \
                        help="On-pulse region. Two values should be provided, " \
                            "the starting phase and the ending phase.")
    parser.add_argument('--sefd', dest='sefd', type=float, \
                        help="The SEFD of the observing system.")
    args = parser.parse_args()
    main()
