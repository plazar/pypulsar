#!/usr/bin/env python
# Simple script to find minimum companion mass for a binary pulsar
# given the mass function
#           Patrick Lazarus, Sept. 18th, 2009

import sys
import optparse

import numpy as np

parser = optparse.OptionParser()
parser.add_option('-m', '--pulsar-mass', dest='mp', type='float', \
                    help="Use given pulsar mass, measured in solar masses." \
                         "(Default: 1.4 solar masses.)", default=1.4)
parser.add_option('-f', '--mass-function', dest='mf', type='float', \
                    help="Use given mass function, measured in solar masses.")
parser.add_option('-i', '--inclination', dest='inclination', type='float', \
                    help="Use given inclination angle (in degrees. Default: 90)", default=90)
(options, sys.argv) = parser.parse_args()

if not hasattr(options, 'mf'):
    sys.stderr.write("Mass function (in solar masses) is required. Exiting!\n")
    sys.exit(1)
mf = options.mf
# Default pulsar mass
mp = options.mp
# Default inclination
sini = np.sin(options.inclination*np.pi/180.0)
if options.inclination > 90 or options.inclination < 0:
    raise ValueError("Inclination angle must be between 0 and 90.")

# Coefficients for cubic equation
a = -mf/sini**3
b = -2*mf*mp/sini**3
c = -mf*mp*mp/sini**3
p = np.array([1,a,b,c])
roots = np.roots(p)
realroots = np.real(roots[np.isreal(roots)])
if realroots.size == 1:
    print "Minimum companion mass (assuming Mp=%g, i=%g): %f Msun" % \
                (options.mp, options.inclination, realroots[0])
else:
    print "Minimum companion mass (assuming Mp=%g, i=%g): " % \
                (options.mp, options.inclination)
    print "\t** Multiple real-valued solutions **"
    print "\t%f Msun" % realroots
    
