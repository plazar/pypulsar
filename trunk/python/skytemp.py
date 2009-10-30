#!/usr/bin/env python
"""
skytemp.py

Program/Module to get sky temperature at different positions
and at different radio frequencies (assuming synchrotron 
spectrum.)

Uses HEALPix Haslam et al. (1982) map from LAMBDA.
(http://lambda.gsfc.nasa.gov/product/forground/fg_haslam_get.cfm)

Requires healpy, which requires pyfits, numpy, and matplotlib:
http://code.google.com/p/healpy/
(pyfits: http://www.stsci.edu/resources/software_hardware/pyfits)
"""

import sys

import healpy
import numpy as np
import matplotlib.pyplot as plt


HASLAM_MAP_FILENAME = "/homes/borgii/alfa/research/PALFA/skytemp/" \
                        "lambda_haslam408_dsds.fits"

HASLAM_FREQ = 408.0 # MHz
SYNCHROTRON_INDEX = -2.7

# Read Haslam map
# healpy outputs some info to stdout, so redirect to /dev/null
stdout = sys.stdout
sys.stdout = open('/dev/null', 'w')
HASLAM_MAP = healpy.read_map(HASLAM_MAP_FILENAME)
sys.stdout = stdout

DEGTORAD = np.pi/180.0

def get_skytemp(gal_long, gal_lat, freq=HASLAM_FREQ, index=SYNCHROTRON_INDEX):
    """Given galactic longitude, latitude and observing
        frequency return the sky temperature in Kelvin
        using the Haslam et al. map.

        Notes: * galactic longitude and latitude should be 
                  provided in degrees.
               * frequency should be provided in MHz.
        
        Optional argument:
        index: synchrotron spectrum index
    """
    
    # Theta and Phi angles defined as expected by healpy
    # Theta and Phi are given expected by healpy in radians
    # In (theta, phi) the centre of Galaxy is (pi/2, 0)
    theta = (90.0 - gal_lat)*DEGTORAD
    phi = (gal_long)*DEGTORAD

    # Read map
    temp_408 = healpy.get_interp_val(HASLAM_MAP, theta, phi)

    # Adjust temperatures for given observing frequency
    temp_obs = change_obsfreq(temp=temp_408, oldfreq=HASLAM_FREQ, \
                                newfreq=freq, index=index)

    return temp_obs


def show_temp_map(freq=HASLAM_FREQ, index=SYNCHROTRON_INDEX, min=None, max=None, \
                    galcoords=[]):
    """Show mollweide projection of Haslam map at frequency 'freq',
        using power law index 'index'. Display temperatures in range
        ('min', 'max'). Temperatures outside this range will be clipped
        to either 'min', or 'max'.
        Mark coordinates provided in 'galcoords'. Entries in list should
        be tuples = (l,b). l and b should be given in degrees.
    """
    healpy.mollview(change_obsfreq(HASLAM_MAP, HASLAM_FREQ, freq, index), \
                        min=min, max=max, unit='K', \
                        title="All-sky map at %gMHz. " \
                                "Haslam et al. (1982) from LAMBDA." % freq)
    if galcoords:
        l,b = np.array(galcoords).transpose()
        healpy.projscatter(l, b, lonlat=True, marker='x', s=50, lw=1.5)
    def keypress(event):
        if event.key in ('q', 'Q'):
            plt.close()
            
    plt.gcf().canvas.mpl_connect('key_press_event', keypress)
    plt.show()


def change_obsfreq(temp, oldfreq, newfreq, index=SYNCHROTRON_INDEX):
    """Change temperature assuming a power law index
        of 'index', using 'oldfreq' and 'newfreq'.
    """
    return temp*(newfreq/oldfreq)**SYNCHROTRON_INDEX


def main():
    deg_symb = unichr(176).encode("latin-1")
    
    galcoords = []
    if options.galcoords:
        for l,b in [g.split(',') for g in options.galcoords]:
            galcoords.append((float(l), float(b)))
        for l, b in galcoords:
            temp = get_skytemp(l, b, freq=options.freq, index=options.index)
            print "l=%g%s, b=%g%s, freq=%g, (index=%g): temp=%g K" % \
                (l, deg_symb, b, deg_symb, options.freq, options.index, temp)
    else:
        sys.stderr.write("No coords provided on command line!\n")
    
    if options.show:
        show_temp_map(freq=options.freq, index=options.index, \
                        min=options.min, max=options.max, \
                        galcoords=galcoords)

if __name__=='__main__':
    import optparse
    parser = optparse.OptionParser(usage="%prog -g coords [-g coords ...] [options]", \
                                    description="Compute and return sky temperature " \
                                                "at given coordiates using Haslam " \
                                                "et al. all-sky map at 408 MHz.", \
                                    version="%prog v0.9 (by Patrick Lazarus")
    parser.add_option('-f', '--freq', dest='freq', type='float', \
                        help="Observing frequency in MHz. (Default: %g MHz.)" % \
                        HASLAM_FREQ, default=HASLAM_FREQ)
    parser.add_option('-i', '--index', dest='index', type='float', \
                        help="Spectral index to use for power law when " \
                             "using a different observing frequency. " \
                             "(Default: %g.)" % SYNCHROTRON_INDEX, \
                             default=SYNCHROTRON_INDEX)
    parser.add_option('-g', '--galcoords', dest='galcoords', type='string', \
                        action='append', \
                        help="Galactic coords to return temperature for. " \
                             "Values should be given in degrees. " \
                             "Longitude and latitude should be seperated " \
                             "by a comma (no space). Multiple -g/--galcoords " \
                             "options are permitted, at least 1 is required.")
    parser.add_option('-s', '--show-map', dest='show', action='store_true', \
                        help='Show map.', default=False)
    parser.add_option('--min', dest='min', type='float', \
                        help="Min temp (in Kelvin) for display range when " \
                             "showing map. All values smaller are shown " \
                             "as same colour as 'min'.", default=None)
    parser.add_option('--max', dest='max', type='float', \
                        help="Max temp (in Kelvin) for display range when " \
                             "showing map. All values larger are shown " \
                             "as same colour as 'max'.", default=None)
    options, args = parser.parse_args()
    main()
