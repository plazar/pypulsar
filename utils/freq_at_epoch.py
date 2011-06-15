#!/usr/bin/env python

"""
Given a parfile and an epoch compute the frequency and
its uncertainty at a requested epoch.
"""
import sys
import numpy as np
import parfile
import psr_utils

print sys.argv

par = parfile.psr_par(sys.argv[1])

for epoch in sys.argv[2:]:
    epoch = float(epoch)
    interval_seconds = (epoch - par.PEPOCH)*psr_utils.SECPERDAY 
    f = par.F0 + interval_seconds*par.F1
    ferr = np.sqrt(par.F0_ERR**2+interval_seconds**2*par.F1_ERR**2)
    print "MJD: %f\n\tf: %0.10f\n\t+- %0.12f" % (epoch, f, ferr)

