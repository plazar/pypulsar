#!/usr/bin/env python

# A program to find a pulsar's position from
# gridding observation
#
#   Fit beam profile to snrs provided (via .pfd files)
#   (See http://www.scipy.org/Cookbook/FittingData "3.2 Fitting a 2D Gaussian")
#
#       Patrick Lazarus, March 18th, 2009

import sys
import numpy as np
import scipy.optimize as opt
from scipy.integrate import trapz
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import prepfold
import psr_utils
from pypulsar.utils.astro import sextant
from pypulsar.utils.astro import protractor
from pypulsar.utils import estimate_snr
debug = 1

class Observation:
    """ Observation object
    """
    def __init__(self, pfdfn):
        """
        Read pfd file and return arrays of snr, ra and dec
    
        Input: pfdfn - a pfd filename
        """
        self.p = prepfold.pfd(pfdfn)
        # RA in arcmin
        print self.p.rastr
        self.ra = protractor.convert(self.p.rastr, 'hmsstr', 'deg')*60
        # Dec in arcmin
        print self.p.decstr
        self.dec = protractor.convert(self.p.decstr, 'dmsstr', 'deg')*60
        self.snr = None
        self.width = None
        
        self.get_snr(self.p)

    def get_snr(self, pfd):
        """
        Return signal-to-noise ratio for the pfd object.

        Input: pfd - a pfd object.
        """
        self.fig = plt.figure()
        self.ax = plt.gca()
        plt.plot(pfd.bestprof.profile, 'k-')
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
            if xmax > len(self.p.bestprof.profile)-1:
                xmax = len(self.p.bestprof.profile)-1
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
        if xmax > len(self.p.bestprof.profile)-1:
            xmax = len(self.p.bestprof.profile)-1
        xmin = int(np.round(xmin))
        xmax = int(np.round(xmax))

        # Calculate S/N using eq. 7.1 from Lorimer and Kramer
        offpulse = np.concatenate((self.p.bestprof.profile[:xmin], \
                                    self.p.bestprof.profile[xmax:]))
        onpulse = self.p.bestprof.profile[xmin:xmax]
        mean = offpulse.mean()
        std = offpulse.std()
        self.oldsnr = ((self.p.bestprof.profile - mean) / std).max()
        weq = trapz((onpulse-mean)/std)/np.max((onpulse-mean)/std)
        self.width = xmax-xmin
        self.snr = np.sum(self.p.bestprof.profile-mean)/std/np.sqrt(weq)

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

def observed_snr(psrsnr, psrra, psrdec):
    """
    Given an intrisic snr, RA, dec and the telescopes beam profile
    return a function that given a beam centre will return the snr
    that will be observed.
    
    psrra and psrdec are in radians
    beam_profile is a function that accepts an angular_separation in radians
    and returns the gain of the beam.
    """
    global beam_profile
    return lambda obsra,obsdec: psrsnr*beam_profile.gain_at_angular_offset(angsep_arcmin(psrra, psrdec, obsra, obsdec))
    
def fit(init_params, data):
    """
    returns (psrsnr, psrra, psrdec), the pulsars parameters found by a fit to data provided
    """
    snrs, ras, decs = data.transpose()
    errorfunction = lambda p: np.ravel(observed_snr(*p)(ras, decs) - snrs)
    if debug:
        output = opt.leastsq(errorfunction, init_params, maxfev=10000, full_output=1)
        print output[-1], output[-2]
        p = output[0]
    else:
        p, success = opt.leastsq(errorfunction, init_params, maxfev=10000)
    return p

def angsep_arcmin(ra1, dec1, ra2, dec2):
    """
    angsep_arcmin(ra1, dec1, ra2, dec2):

    return angular separation in arcmin.
    ra1, dec1, ra2, dec2 are all given in armin.
    """
   
    ra1_rad = ra1/60.0*psr_utils.DEGTORAD
    ra2_rad = ra2/60.0*psr_utils.DEGTORAD
    dec1_rad = dec1/60.0*psr_utils.DEGTORAD
    dec2_rad = dec2/60.0*psr_utils.DEGTORAD
    angsep_rad = np.arccos(np.sin(dec1_rad)*np.sin(dec2_rad)+np.cos(dec1_rad)*np.cos(dec2_rad)*np.cos(ra1_rad-ra2_rad))
    angsep_arcmin = angsep_rad*psr_utils.RADTODEG*60.0
    return angsep_arcmin


def main():
    observations = [Observation(pfdfn) for pfdfn in sys.argv[1:]]
    data = np.array([(o.snr, o.ra, o.dec) for o in observations])
    if debug:
        print "data:"
        for z in data:
            print "\tSNR:", z[0], "RA:", z[1], "Dec:", z[2]
    # Init guess is max SNR, weighted avg of RA, weighted avg of Dec
    data_T = data.transpose()
    init_params = (data_T[0].max(), (data_T[0]*data_T[1]).sum()/data_T[0].sum(), \
                        (data_T[0]*data_T[2]).sum()/data_T[0].sum())
    if debug:
        print "initial parameters:"
        print "\tSNR:", init_params[0], "RA:", init_params[1], "Dec:", init_params[2]
    global beam_profile
    # Use gain = 1
    beam_profile = estimate_snr.EstimateFWHMSNR(3.35/2.0, 1420, 100, 2, 1, 24)
    result = fit(init_params, data) 
    if debug:
        print "results:"
        print "\tSNR:", result[0], "RA:", result[1], "Dec:", result[2]
    psrsnr, psrra, psrdec = result
    snrs, ras, decs = data.transpose()
    plt.figure(figsize=(8.5,11))
    plt.subplot(211)
    plt.title("Fitting gridding observations to determine pulsar position")
    plt.scatter((ras-psrra)*60/15.0, (decs-psrdec)*60, c=snrs, marker='o', label='_nolegend_')
    plt.spring()
    cbar = plt.colorbar()
    cbar.set_label(r"$SNR$")
    plt.scatter(np.array([0]), np.array([0]), s=100, c='k', marker=(5,1,0), \
                    label='Best PSR posn')
    if debug:
        plt.scatter(np.array([init_params[1]-psrra])*60/15.0, \
                        np.array([init_params[2]-psrdec])*60, \
                            s=100, c='w', marker=(5,1,0), label='Init PSR posn')
    plt.legend(loc='best')
    plt.xlabel("RA (sec) + %02.0f:%02.0f:%07.4f" % psr_utils.rad_to_hms(psrra/60.0*psr_utils.DEGTORAD))
    plt.ylabel("Dec (arcsec) + %02.0f:%02.0f:%07.4f" % psr_utils.rad_to_dms(psrdec/60.0*psr_utils.DEGTORAD))
   
    obsangseps = np.zeros(len(snrs))
    for ii in range(len(snrs)):
        obsangseps[ii] = angsep_arcmin(psrra, psrdec, ras[ii], decs[ii])
    maxangsep = obsangseps.max()
    angseps = np.linspace(0,maxangsep*1.1, 1000)
    plt.subplot(212)
    plt.plot(angseps, psrsnr*beam_profile.gain_at_angular_offset(angseps), 'k', zorder=-1)
    plt.scatter(obsangseps, snrs, c=snrs, zorder=1)
    plt.xlabel("Angular separation (arcmin)")
    plt.ylabel("SNR")
    plt.savefig('gridding.tmp.ps', papertype='letter', orientation='portrait')
    cid_keypress = plt.gcf().canvas.mpl_connect('key_press_event', \
                                                    keypress)
    plt.show()


def keypress(event):
    if event.key in ('q', 'Q'):
        sys.exit(0)


if __name__ == '__main__':
    main()
