# simple script to solve mass fuction assuming sini=1, m_pulsar=1.4
import sys
import numpy as np

if len(sys.argv) <= 1:
    sys.stderr.write("Need to provide mass function (in solar masses) as first argument. Exiting!")
    sys.exit(1)
mf = float(sys.argv[1])
# Default pulsar mass
mp = 1.4
# Default inclination
sini = 1

# Coefficients for cubic equation
a = -mf/sini**3
b = -2*mf*mp/sini**3
c = -mf*mp*mp/sini**3
p = np.array([1,a,b,c])
roots = np.roots(p)
realroots = np.real(roots[np.isreal(roots)])
if realroots.size == 1:
    print "Minimum pulsar mass (assuming Mp=1.4):", realroots[0], "M_solar"
else:
    print "Minimum pulsar mass (assuming Mp=1.4): ** Multiple real-valued solutions **"
    print realroots, "M_solar"
    
