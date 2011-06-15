#!/usr/bin/env python
import glob
import os.path

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

import presto
import infodata

FUDGEFACTOR = 1

class FreqInterval:
    def __init__(self, fcent, ferr, numel=1):
	self.ferr = ferr
	self.fcent = fcent
	self.flo = self.fcent - self.ferr
	self.fhi = self.fcent + self.ferr
	self.width = (self.fhi-self.flo)*FUDGEFACTOR
	self.numelements = numel

    def __contains__(self, other):
	if type(self) != type(other):
	   raise ValueError("Contains test must be made between two FreqInterval objects.")
	return ((self.flo) < other.flo < (self.fhi)) or \
		((self.flo) < other.fhi < (self.fhi)) or \
		((other.flo) < self.flo < (other.fhi)) or \
		((other.flo) < self.fhi < (other.fhi))

    def __add__(self, other):
	if type(self) != type(other):
	    raise ValueError("Addition must be between two FreqInterval objects.")
	flo = min([self.fcent-self.ferr, other.fcent-other.ferr])
	fhi = max([self.fcent+self.ferr, other.fcent+other.ferr])
	return FreqInterval((flo+fhi)/2.0, (fhi-flo)/2.0, \
			    numel=self.numelements+other.numelements)

    def __str__(self):
	return "<FreqInterval: flo=%g, fhi=%g, numelements=%d>" % \
		(self.flo, self.fhi, self.numelements)

    def zaplist_string(self):
	return "\t%f\t%f" % (self.fcent, self.width)

def main():
    freqs = []
    freqerrs = []
    filenums = []
    intervals = []

    filenum = 0
    inffiles = glob.glob("*.inf")
    for inffile in inffiles:
        accelfile = inffile[:-4]+"_ACCEL_0.cand"
   	if not os.path.exists(accelfile):
	    continue
	filenum += 1
    	rzws = presto.read_rzwcands(accelfile)
    	inf = infodata.infodata(inffile)
    	T = inf.dt*inf.N
    	for rzw in rzws:
	    freq = rzw.r/T
	    freqerr = rzw.rerr/T
	    freqs.append(freq)
	    freqerrs.append(freqerr)
	    filenums.append(filenum)

	    fint = FreqInterval(freq, freqerr)
	    # Traverse list of intervals backwards
	    for ii in range(len(intervals))[::-1]:
		if fint in intervals[ii]:
		    fint = fint+intervals[ii]
		    matchfound = True
		    del intervals[ii]
	    intervals.append(fint)

    freqs = np.array(freqs)
    freqerrs = np.array(freqerrs)
    filenums = np.array(filenums)

    plt.figure(figsize=(11,8.5))
    ebax = plt.axes((0.1, 0.1, 0.7, 0.7))
    plt.errorbar(freqs, filenums, xerr=freqerrs, fmt=None, zorder=1, ecolor='k')
    # Plot intervals worth masking
    for i in intervals:
	if i.numelements > 7:
	    r = matplotlib.patches.Rectangle((i.fcent-i.width/2.0,0), i.width, max(filenums), \
						fill=True, fc='r', ec='none', \
						alpha=0.25, zorder=-1)
	    plt.gca().add_patch(r)
	    print i.zaplist_string()
    plt.xlabel("Spin Frequency (Hz)")
    plt.ylabel("File number (index)")
    hax = plt.axes((0.8, 0.1, 0.15, 0.7), sharey=ebax)
    plt.hist(filenums, bins=max(filenums), range=(0, max(filenums)), orientation='horizontal', fc='none')
    plt.savefig("accelcands.ps", orientation="landscape", papertype="letter")
    plt.show()

if __name__=='__main__':
    main()
