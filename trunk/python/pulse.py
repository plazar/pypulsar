#!/usr/bin/env python

"""
pulse.py

Pulsar pulse object.

Patrick Lazarus, June 24, 2009
"""
import copy
import numpy as np
import scipy.signal

class Pulse:
    def __init__(self, number, mjd, time, durration, profile, on_pulse_regions=None):
        """Create a pulse object. Arguments provided are:
            - number: The pulse number (for example, counted from
                        beginning of an observation.)
            - mjd: MJD of the first sample in the pulse profile.
            - time: Time elapsed (in seconds) measured from
                        beginning of observation to first sample
                        in the pulse profile.
            - duration: Durration of the pulse (in seconds).
            - profile: A numpy array containing the raw pulse
                        profile (a slice of a timeseries).
            - on_pulse_regions: A list of 2-tuples. Each tuple
                        defines an on-pulse region measured in
                        rotational phase (between 0.0 and 1.0).
        """
        self.number = number
        self.mjd = mjd
        self.time = time
        self.duration = duration
        self.profile = profile.flatten()
        self.N = profile.size # number of samples in the profile
        self.set_on_pulse_regions(on_pulse_regions) # list of 2-tuples defining
                                                    # on-pulse regions


    def __str__(self):
        return "Pulse #: %d\n\tMJD: %0.15f\n\tTime: %8.2f\n\tDuration: %8.4f" % \
                    (self.number, self.mjd, self.time, self.duration)


    def set_onoff_pulse_regions(self, on_pulse_regions):
        """Set the on-pulse and off-pulse regions using
            on-pulse regions provided.
            'on_pulse_regions': is a list of 2-tuples.
            Each 2-tuple contains the start and end of
            a on-pulse region. The values in the tuples
            are measured in fraction of rotational phase
            (between 0.0 and 1.0).
        """
        on_pulse = np.array(on_pulse_regions).astype('float32')
        # Sort based on first value of each pair
        on_pulse = on_pulse[on_pulse.argsort(axis=0).transpose()[0]]
        # Check if there are any overlaps (or inverted pairs)
        if np.any(on_pulse.flat != sorted(on_pulse.flat)):
            raise OnPulseRegionError("On-pulse regions overlap or are inverted")
        else:
            # Everything checks out
            self.on_pulse = on_pulse
            
            # Set off-pulse regions.
            # Be careful if on-pulse starts or ends at beginning
            # or end of rotational phase.
            off_pulse = self.on_pulse.flatten()
            if off_pulse[0] == 0.0:
                off_pulse = off_pulse[1:]
            else:
                off_pulse = np.concatenate(([None], off_pulse))
            if flattened_on_pulse[-1] == 1.0:
                off_pulse = off_pulse[:-1]
            else:
                off_pulse = no.concatenate((off_pulse, [None]))
            off_pulse.shape = (off_pulse.size/2, 2)
            self.off_pulse = off_pulse


    def get_data(self, regions):
        """Return numpy array of data from 'regions' concatenated
            together. 'regions' is a list of 2-tuples. Each 2-tuples
            contains a region of data to extract measured in
            rotational phase (between 0.0 and 1.0).
        """
        data = []
        for (lo, hi) in regions:
            lobin = int(self.N * lo)
            hibin = int(self.N * hi)
            if hibin >= lobin:
                raise OnPulseRegionError("lobin=%s, hibin=%s" % (lobin, hibin))
            data.append(self.profile[lobin:hibin])
        return np.concatenate(data)


    def get_on_pulse(self):
        """Return numpy array of on-pulse regions concatenated
            together.
        """
        return get_data(self.on_pulse)
    
    def get_off_pulse(self):
        """Return numpy array of off-pulse regions concatenated
            together.
        """
       return get_data(self.off_pulse) 
        off_pulse_data = [


    def make_copy(self):
        """Return a deep copy of 'self'
            (Uses Python's copy.deepcopy(...)).
        """
        return copy.deepcopy(self)


    def scale(self):
        """Scale 'self'.
            Scaling subtracts the mean of the off-pulse regions
            and divides by the stddev of the off-pulse regions.

            NOTE: The profile attribute of 'self' will be modified.
        """
        off_pulse_region = self.get_off_pulse()
        # Do scaling
        self.profile -= np.mean(off_pulse_region)
        self.profile /= np.std(off_pulse_region)
    
    
    def downsample(self, downfactor=1):
        """Downsample profile by adding 'downfactor' adjacent
            bins. ('downfactor' does not have to be a factor
            of the size of the profile. in this case the end
            of the profile will be truncated.)

            NOTE: The profile attribute of 'self' will be modified.
        """
        if downfactor > 1:
            self.N = self.N/downfactor*downfactor # New length of profile
            self.profile = self.profile[:self.N]
            self.profile.shape = (self.N/downfactor, downfactor)
            self.profile = self.profile.sum(axis=1)
    

    def smooth(self, smoothfactor=1):
        """Smooth profile by convolving with tophat of width
            'smoothfactor'. The height of the tophat is chosen
            such that RMS = 1 after smoothing. No overlap is
            used so the edges of the profile might have less
            signal than expected.

            This bit of code is taken from Scott Ransom's
            PRESTO's single_pulse_search.py (line ~ 423).

            NOTE: The profile attribute of 'self' will be modified.
        """
        if smoothfactor > 1:
            kernel = np.ones(smoothfactor, dtype='float32') / \
                        np.sqrt(smoothfactor)
            #
            # Do we want to wrap profile around to avoid edge
            # effects of the convolution?
            #
            self.profile = scipy.signal.convolve(self.profile, kernel, 'same')


    def detrend(self, numchunks=5):
        """Break profile into 'numchunks' and remove a linear trend
            from each chunk.

            NOTE: The profile attribute of 'self' will be modified.
        """
        break_points = np.round(np.linspace(0, self.N, numchunks+1))
        self.profile = scipy.signal.detrend(self.profile, bp=break_points)
        

    def is_masked(self, numchunks=5):
    """Break pulse profile into 'numchunks'. Check each chunk
        if it is flat (max value == min value). If it is then
        profile is partially masked, return True. Otherwise
        return False.
    """
    edges = np.round(np.linspace(0, self.profile.size, numchunks+1))
    for i in range(0,N):
        if self.profile[edges[i]:edges[i+1]].ptp() == 0:
            # Current section of profile is flat.
            # Profile is partially flat.
            return True
    # No section of profile is flat.
    return False


class OnPulseRegionError(Exception):
    """Error when on-pulse region is ill-defined.
    """
    def __init__(self, message):
        self.message = message


    def __str__(self):
        return "On-pulse region is ill-defined.", message
