#!/usr/bin/env python

import sys

from pypulsar.utils.astro import sextant

print sextant.equatorial_to_galactic(float(sys.argv[1]), float(sys.argv[2]), input='deg')
