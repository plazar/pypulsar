import numpy as np

"""
Zenith angle dependence of the L-wide receiver at Arecibo.

Data are taken from http://www2.naic.edu/~phil/lbw/gainMar03/lbwGainCurvesMar03.html
Specifically read from plot http://www.naic.edu/~phil/lbw/gainMar03/lbwgainfitMar03.ps (Freq=1550 MHz)

Patrick Lazarus Jan 15, 2014
"""

def gain(za, az):
    za_less_14 = np.clip(za-14, 0, np.inf)
    val = 10.14891 + 0.03814*za - 0.05113 * za_less_14**2 - 0.00193*za_less_14**3
    return val


def tsys(za, az):
    return 30.0
