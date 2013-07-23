#!/usr/bin/env python
"""
Manually determine on-pulse region for the profile in a pfd file
    and calculate the SNR.

    Patrick Lazarus, May 14, 2013 (taken out of gridding.py).
"""

import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.integrate import trapz

import prepfold
import psr_utils

from pypulsar.utils.astro import protractor
debug = 1

class Observation:
    """ Observation object
    """
    def __init__(self, pfdfn):
        """
        Read pfd file and return arrays of snr, ra and dec
    
        Input: pfdfn - a pfd filename
        """
        self.fn = pfdfn
        self.p = prepfold.pfd(pfdfn)
        # RA in arcmin
        self.ra = protractor.convert(self.p.rastr, 'hmsstr', 'deg')*60
        # Dec in arcmin
        self.dec = protractor.convert(self.p.decstr, 'dmsstr', 'deg')*60
        self.snr = None
        self.width = None
        
        prof = self.p.bestprof.profile
        imax = np.argmax(prof)
        nrot = (imax-len(prof)/2) % len(prof)
        print "Profile maximum at bin %d. Rotating by %d bins." % (imax, nrot)
        self.prof = psr_utils.rotate(prof, nrot)
        
        self.get_snr()

    def get_snr(self):
        """
        Return signal-to-noise ratio for the pfd object.

        Input: pfd - a pfd object.
        """
        self.fig = plt.figure()
        self.ax = plt.gca()
        plt.plot(self.prof, 'k-')
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

        # Calculate S/N using eq. 7.1 from Lorimer and Kramer
        offpulse = np.concatenate((self.prof[:xmin], \
                                    self.prof[xmax:]))
        onpulse = self.prof[xmin:xmax]
        mean = offpulse.mean()
        std = offpulse.std()
        self.oldsnr = ((self.prof - mean) / std).max()
        weq = trapz((onpulse-mean)/std)/(np.max((onpulse-mean)/std))
        self.width = xmax-xmin
        self.snr = np.sum(self.prof-mean)/std/np.sqrt(weq)
        if debug:
            print "Width selected (bins):", self.width
            print "Equivalent width (bins):", weq
            print "Integral under the normalised on-pulse region:", \
                    trapz((onpulse-mean)/std)
            print "SNR:", self.snr

        # disconnect mousemove, mouserelease and reset state
        self.eventpress = None
        self.mousepressed = False
        self.fig.canvas.mpl_disconnect(self.cid_mousemove)
        self.fig.canvas.mpl_disconnect(self.cid_mouserelease)

        # Display width and snr
        print "SNR:", self.snr
        print "Width (bins):", self.width

    def keypress(self, event):
        if event.key == 'enter' and self.snr is not None and \
                self.width is not None:
            plt.close(event.canvas.figure.number)


def main():
    pfdfns = sys.argv[1:]
    for pfdfn in pfdfns:
        print pfdfn
        obs = Observation(pfdfn)        
        print "%s \tSNR: %s" % (obs.fn, obs.snr)


if __name__ == '__main__':
    main()
