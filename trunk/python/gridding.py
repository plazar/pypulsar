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
import pylab as plt
import prepfold
import psr_utils
import coordconv
import estimate_snr
debug = 1

def read_pfds(pfdfns):
    """
    Read pfd files and return arrays of snr, ra and dec
    
    Input: pfdfns - a list of pfd filenames
    """
    N = len(pfdfns)
    snrs = np.zeros(N)
    ras = np.zeros(N)
    decs = np.zeros(N)
    for ii, pfdfn in enumerate(pfdfns):
        p = prepfold.pfd(pfdfn)
        snrs[ii] = ((np.array(p.bestprof.profile)-p.bestprof.prof_avg)/p.bestprof.prof_std).max()/np.sqrt(p.T/p.bary_p1)
        ras[ii] = coordconv.rastr_to_deg(coordconv.fmrastr_to_rastr(p.rastr))*60 # RA in arcmin
        decs[ii] = coordconv.decstr_to_deg(coordconv.fmdecstr_to_decstr(p.decstr))*60 # dec in arcmin
    return snrs, ras, decs

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
    snrs, ras, decs = data
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
    data = read_pfds(sys.argv[1:])
    snrs, ras, decs = data
    if debug:
        print "data:"
        for z in zip(snrs, ras, decs):
            print "\tSNR:", z[0], "RA:", z[1], "Dec:", z[2]
    init_params = (snrs.max(), ras.mean(), decs.mean())
    if debug:
        print "initial parameters:"
        print "\tSNR:", init_params[0], "RA:", init_params[1], "Dec:", init_params[2]
    global beam_profile
    beam_profile = estimate_snr.EstimateFWHMSNR(3.35/2.0, 1420, 100, 2, 10.4, 24)
    result = fit(init_params, data) 
    if debug:
        print "results:"
        print "\tSNR:", result[0], "RA:", result[1], "Dec:", result[2]
    psrsnr, psrra, psrdec = result
    plt.figure(figsize=(8.5,11))
    plt.subplot(211)
    plt.title("Fitting gridding observations to determine pulsar position")
    plt.scatter((ras-psrra)*60/15.0, (decs-psrdec)*60, c=snrs/psrsnr, marker='o', label="_nolegend_")
    plt.scatter(np.array([0]), np.array([0]), c='1', s=36, marker='d', label="Pulsar Position")
    plt.xlabel("RA (sec) + %02.0f:%02.0f:%07.4f" % psr_utils.rad_to_hms(psrra/60.0*psr_utils.DEGTORAD))
    plt.ylabel("Dec (arcsec) + %02.0f:%02.0f:%07.4f" % psr_utils.rad_to_dms(psrdec/60.0*psr_utils.DEGTORAD))
    plt.legend(loc='best')
   
    obsangseps = np.zeros(len(snrs))
    for ii in range(len(snrs)):
        obsangseps[ii] = angsep_arcmin(psrra, psrdec, ras[ii], decs[ii])
    maxangsep = obsangseps.max()
    angseps = np.linspace(0,maxangsep*1.1, 1000)
    plt.subplot(212)
    plt.plot(angseps, psrsnr*beam_profile.gain_at_angular_offset(angseps), 'k')
    plt.scatter(obsangseps, snrs, c='k')
    plt.xlabel("Angular separation (arcmin)")
    plt.ylabel("SNR")
    plt.savefig('gridding.tmp.ps', papertype='letter', orientation='portrait')
    plt.show()

if __name__ == '__main__':
    main()
