"""
datfile.py

PRESTO .dat file object.

Patrick Lazarus, June 16, 2009
"""

import os.path
import numpy as np
import psr_utils
import infodata
import pulse

# Define constants and default values
DTYPE = 'float32'       # binary data in PRESTO's .dat files 
                        #   is stored as 32-bit floats

class Datfile:
    def __init__(self, datfn, dtype=DTYPE):
        self.datfn = datfn
        self.dtype = np.dtype(dtype)
        if self.datfn.endswith(".dat"):
            self.basefn = self.datfn[:-4]
            self.datfile = open(datfn, 'r')
            self.inffn = "%s.inf" % self.basefn
            self.infdata = infodata.infodata(self.inffn)
        else:
           raise "Filename (%s) doesn't end with '.dat'"
        # initialize file, and current sample and current time counters
        self.rewind()

        
    def __str__(self):
        return self.datfn


    def read_Nsamples(self, N):
        """Read N samples from datfile and return a numpy array
            of the data.
        """
        new_curr_sample = self.currsample + N
        if new_curr_sample > self.infdata.N:
            return None
        else:
            self.currsample = new_curr_sample
            return np.fromfile(self.datfile, dtype=self.dtype, count=N)


    def read_Tseconds(self, T):
        """Read T seconds worth of data from datfile and return
            a numpy array of the data.
        """
        # Compute number of samples to read
        endsample = np.round((self.currtime+T)/self.infdata.dt)
        num_samples_to_read = int(endsample-self.currsample)
        # Update current time
        self.currtime += T
        # Read data using read_Nsamples(...)
        return self.read_Nsamples(num_samples_to_read)


    def rewind(self):
        """Rewind file to beginning. Also, reset current time
            and current sample.
        """
        self.datfile.seek(0)
        self.currsample = 0
        self.currtime = 0.0


    def pulses(self, period_at_mjd):
        """Generator method to generate/iterate over pulse
            profiles from datfile.
            'period_at_mjd' is a function that, given a mjd,
            will return the spin period at that mjd.
        """
        # Initialize
        self.rewind()
        pulse_number = 1
        current_mjd = self.infdata.epoch    # MJD at start of obs
        current_obstime = 0                 # Time into obs (seconds)
        current_period = period_at_mjd(current_mjd)
        current_pulse = self.read_Tseconds(current_period)
        while current_pulse is not None:
            # yield (return) current pulse (and other information)
            yield pulse.Pulse(pulse_number, current_mjd, current_obstime, \
                    current_period, current_pulse)
            # Update pulse number, mjd, period and pulse
            pulse_number += 1
            current_mjd += float(current_period) / psr_utils.SECPERDAY
            current_obstime += current_period
            current_period = period_at_mjd(current_mjd)
            current_pulse = self.read_Tseconds(current_period)
