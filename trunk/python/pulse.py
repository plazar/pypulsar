#!/usr/bin/env python

"""
pulse.py

Pulsar pulse object.

Patrick Lazarus, June 24, 2009
"""
import copy
import os.path
import types
import numpy as np
import scipy.signal

# Import matplotlib/pyplot and set for non-interactive plots
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

class Pulse:
    def __init__(self, number, mjd, time, duration, profile, \
                    origfn, dt, on_pulse_regions=None):
        """Create a pulse object. Arguments provided are:
            - number: The pulse number (Counted from beginning 
                                        of an observation.)
            - mjd: MJD of the first sample in the pulse profile.
            - time: Time elapsed (in seconds) measured from
                        beginning of observation to first sample
                        in the pulse profile.
            - duration: Duration of the pulse (in seconds).
            - profile: A numpy array containing the raw pulse
                        profile (a slice of a timeseries).
            - origfn: Name of data file the pulse originated from.
            - dt: width of each profile bin (in seconds).
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
        self.dt = dt # width of bin (in seconds)
        self.origfn = origfn
        if (type(on_pulse_regions)==types.ListType or \
                type(on_pulse_regions)==np.ndarray) and len(on_pulse_regions):
            self.set_onoff_pulse_regions(on_pulse_regions)
        else:
            self.on_pulse = None
            self.off_pulse = None

    def __str__(self):
        return "Pulse #: %s\n\tMJD: %0.15f\n\tTime: %8.2f s\n\tDuration: %8.4f s\n" % \
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

            # Flattened version of on_pulse will become off_pulse
            off_pulse = self.on_pulse.flatten() 
            if off_pulse[0] == 0.0:
                off_pulse = off_pulse[1:]
            else:
                off_pulse = np.concatenate(([None], off_pulse))
            if off_pulse[-1] == 1.0:
                off_pulse = off_pulse[:-1]
            else:
                off_pulse = np.concatenate((off_pulse, [None]))
            off_pulse.shape = (off_pulse.size/2, 2)
            self.off_pulse = off_pulse

    def get_data(self, regions=None):
        """Return numpy array of data from 'regions' concatenated
            together. 'regions' is a list of 2-tuples. Each 2-tuples
            contains a region of data to extract measured in
            rotational phase (between 0.0 and 1.0).
        """
        if regions is None or len(regions) == 0:
            regions = [(None, None)] # return all data
        data = []
        for (lo, hi) in regions:
            if lo is None:
                lobin = None
            else:
                lobin = int(self.N * lo)
            if hi is None:
                hibin = None
            else:
                hibin = int(self.N * hi)
            if lobin is not None and hibin is not None and hibin <= lobin:
                raise OnPulseRegionError("lobin=%s, hibin=%s" % (lobin, hibin))
            data.append(self.profile[lobin:hibin])
        return np.concatenate(data)

    def get_on_pulse(self):
        """Return numpy array of on-pulse regions concatenated
            together.
        """
        return self.get_data(self.on_pulse)
    
    def get_off_pulse(self):
        """Return numpy array of off-pulse regions concatenated
            together.
        """
        return self.get_data(self.off_pulse) 

    def get_pulse_energies(self):
        """Return on- and off- pulse energies based on scaled 
            pulse profile. Return value is a tuple:
            (on-energy, off-energy)

            on-pulse energy = sum(on-pulse region)
            off-pulse energy = sum(off-pulse energy)
        """
        copy = self.make_copy()
        copy.scale()
        on_energy = np.sum(copy.get_on_pulse())
        off_energy = np.sum(copy.get_off_pulse())
        return (on_energy, off_energy)

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
            bins. ('downfactor' must be a factor of the size 
            of the profile. Otherwise an exception is thrown.)

            NOTE: The profile attribute of 'self' will be modified.
        """
        if downfactor > 1:
            # downfactor does not need to be a factor of the size
            # of the profile. Left over samples will be truncated.
            # This is not ideal so raise an exception.
            if self.N % downfactor != 0:
                raise "downfactor (%d) is not a factor of profile length (%d)! ... need proper exception." % (downfactor, self.N)
            self.profile = self.profile[:self.N/downfactor*downfactor]
            self.N = int(self.N/downfactor) # New length of profile
            self.profile.shape = (self.N, downfactor)
            self.profile = self.profile.sum(axis=1)
            self.dt *= downfactor 

    def downsample_Nbins(self, N):
        """Downsample profile so final profile has N bins. 

            NOTE: The profile attribute of 'self' will be modified.
        """
        if N > self.N:
            raise "Cannot downsample so new profile (%d) is longer than old profile (%d)! ... need proper exception." % (N, self.N)
        downfactor = int(self.N/N)
        numleftover = self.N % N
        leftover = self.profile[self.N - numleftover:] # should be less than downfactor long
        self.profile = self.profile[:self.N - numleftover]
        self.profile.shape = (N, downfactor)
        self.profile = self.profile.mean(axis=1)
        # DONT Add leftovers to last bin
        # self.profile[-1] += leftover.sum()
        self.N = N # New length of profile
        self.dt *= downfactor
            
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
            # Wrap profile around to avoid edge
            # effects of the convolution
            #
            # NOTE: This isn't ideal because we're not dealing with a
            #       folded profile. Each pulse is a unique piece of the
            #       timeseries.
            prof = np.concatenate([self.profile[-smoothfactor:], \
                                self.profile, self.profile[:smoothfactor]])
            smooth_prof = scipy.signal.convolve(prof, kernel, 'same')
            self.profile = smooth_prof[smoothfactor:-smoothfactor]

    def detrend(self, numchunks=5):
        """Break profile into 'numchunks' and remove a linear trend
            from each chunk.

            NOTE: The profile attribute of 'self' will be modified.
        """
        break_points = np.round(np.linspace(0, self.N, numchunks+1))
        self.profile = scipy.signal.detrend(self.profile, bp=break_points)
        
    def interpolate(self, numsamples):
        """Interpolate profile so it has 'numsamples' across it.

            NOTE: The profile attribute of 'self' will be modified.
        """
        xp = np.arange(self.N)
        x = np.linspace(0, self.N-1, numsamples)
        self.profile = np.interp(x, xp, self.profile)
        self.dt = self.dt*self.N/float(numsamples)
        self.N = numsamples

    def interp_and_downsamp(self, numsamples):
        """Interpolate and then downsample profile so it has
            'numsamples' samples when finished.

            NOTE: The profile attribute of 'self' will be modified.
        """
        downsamp = int(self.N/numsamples)+1
        interp = downsamp*numsamples
        
        # The following settings takes much longer, but is possibly
        # more accurate (when overplotted with above method the two
        # profiles look identical, by eye).
        # downsamp = self.N
        # interp = numsamples * self.N
        
        self.interpolate(interp)
        self.downsample(downsamp)

    def is_masked(self, numchunks=5):
        """Break pulse profile into 'numchunks'. Check each chunk
            if it is flat (max value == min value). If it is then
            profile is partially masked, return True. Otherwise
            return False.
        """
        edges = np.round(np.linspace(0, self.profile.size, numchunks+1))
        for i in range(0,numchunks):
            if self.profile[edges[i]:edges[i+1]].ptp() == 0:
                # Current section of profile is flat.
                # Profile is partially flat.
                return True
        # No section of profile is flat.
        return False

    def plot(self, basefn=None, downfactor=1):
        """Plot the pulse profile.
            'basefn' is base filename to use. Default will
            be to use same as original data file's name
            (with the extension).
            Downsample profiles by factor 'downfactor' before plotting.
        """
        #
        # Assume original filename has an extension
        # NOTE: This is probably not good to assume.
        #
        if basefn is None:
            basefn, extension = os.path.splitext(self.origfn)
        copy_of_self = self.make_copy()
        # Interpolate before downsampling
        interp = ((copy_of_self.N/downfactor)+1)*downfactor
        copy_of_self.interpolate(interp)
        copy_of_self.downsample(downfactor)
        plt.figure()
        plt.plot(copy_of_self.profile, 'k-', lw=0.5)
        plt.xlabel("Profile bin")
        plt.title("Pulse #%d" % self.number)
        plt.savefig("%s.prof%d.ps" % (basefn, self.number), \
                        orientation='landscape')
            

    def write_to_file(self, basefn=None):
        """Dump the pulse to file.
            'basefn' is base filename to use. Default will
            be to use same as original data file's name
            (with the extension).
        """
        #
        # Assume original filename has an extension
        # NOTE: This is probably not good to assume.
        #
        if basefn is None:
            basefn, extension = os.path.splitext(self.origfn)
        file = open("%s.prof%d" % (os.path.split(basefn)[1], self.number), 'w')
        file.write("# Original data file              = %s\n" % self.origfn)
        file.write("# Pulse Number                    = %d\n" % self.number)
        file.write("# MJD of start of pulse           = %0.15f\n" % self.mjd)
        file.write("# Time into observation (seconds) = %f\n" % self.time)
        file.write("# Duration of pulse (seconds)     = %0.15f\n" % self.duration)
        file.write("# Profile bins                    = %d\n" % self.N)
        file.write("# Width of profile bin (seconds)  = %g\n" % self.dt)
        if self.on_pulse is not None:
            for i, (lo,hi) in enumerate(self.on_pulse):
                file.write("# On-pulse region %2d (phase)      = %f-%f\n" % \
                                                                    (i,lo,hi))
        file.write("###################################\n")
        for i, val in enumerate(self.profile):
            file.write("%-10d %f\n" % (i, val))
        file.close()

    def to_summed_pulse(self):
        """Create and return a SummedPulse object using
            self as the starting point of the sum.
        """
        summed_pulse = SummedPulse(self.number, self.mjd, self.time, \
                                self.duration, self.profile, self.origfn, \
                                self.dt, self.on_pulse)
        # summed_pulse.scale() ## DEBUG
        return summed_pulse
        
    def __add__(self, other):
        """Add two Pulse objects.
        """
        if hasattr(other, 'pulse_registry'):
            # Other is a SummedPulse object. Use SummedPulse's __iadd__ method.
            summed_pulse = other.make_copy()
        else:
            copy_of_other = other.make_copy()
            # copy_of_other.scale() ## DEBUG
            summed_pulse = copy_of_other.to_summed_pulse()
        summed_pulse += self
        return summed_pulse


class SummedPulse(Pulse):
    """Object to keep track of sum of multiple pulses.
        Sub-class of Pulse.

        SummedPulse is very similar to Pulse except it
        keeps track of what pulses have been summed.
    """
    def __init__(self, number, mjd, time, duration, profile, \
                    origfn, dt, on_pulse_regions=None, 
                    init_registry=None, init_count=1):
        # Call superclass' constructor
        Pulse.__init__(self, number, mjd, time, duration, profile, \
                    origfn, dt, on_pulse_regions)
        # Initialize registry to keep track of what pulses are summed
        if init_registry is not None:
            self.pulse_registry = init_registry
        else:
            self.pulse_registry = {origfn: [number]}
        self.count = init_count

    def __iadd__(self, other):
        """Sum other with self as self.
            Other can be a SummedPulse or a Pulse.
        """
        if self.dt != other.dt:
            raise "Incompatible binwidths!... need proper exception"
        if hasattr(other, 'pulse_registry'):
            # Merge other's registry into self's registry
            # If there is a conflict, pulse is included in
            # self and other, then raise an Exception.
            for fn in other.pulse_registry.keys():
                if fn in self.pulse_registry.keys():
                    for num in other.pulse_registry[fn]:
                        if num in self.pulse_registry[fn]:
                            # Conflict
                            raise "Conflict exception... need proper exception"
                        else:
                            self.pulse_registry[fn].append(num)
                else:
                    # Use slice to make copy of other's entry in registry
                    self.pulse_registry[fn] = other.pulse_registry[fn][:]
            ocount = other.count
        else:
            # 'other' is a Pulse, not SummedPulse
            if other.origfn in self.pulse_registry.keys():
                self.pulse_registry[other.origfn].append(other.number)
            else:
                self.pulse_registry[other.origfn] = [other.number]
            ocount = 1
       
        # Prepare profiles for summing
        # self.scale()
        copy_of_other = other.make_copy()
        # copy_of_other.scale() ## DEBUG
        
        # Truncate to size of smaller profile
        self.N = np.min([self.N, copy_of_other.N])
        self.duration = np.min([self.duration, copy_of_other.duration])
        self.profile = self.profile[:self.N] + copy_of_other.profile[:self.N]
        
        # Update epoch of pulse
        #
        # NOTE: number is only meaningful if all pulses come from same obs.
        #
        self.number = (self.count * self.number + ocount * other.number) / \
                                float(self.count + ocount)
        self.time = (self.count * self.time + ocount * other.time) / \
                                float(self.count + ocount)
        self.mjd = (self.count * self.mjd + ocount * other.mjd) / \
                                float(self.count + ocount)

        # Update count of profiles summed
        self.count += ocount
        
        #
        # NOTE: Should we modify on-pulse region?
        #
        return self

    def __contains__(self, item):
        """Check if item is contained in self.
            If item is a SummedPulse object return True
            if there is _any_ overlap between registries.
            If item is a Pulse object return True if
            it is listed in self's registry.
        """
        if hasattr(item, 'pulse_registry'):
            for fn in item.pulse_registry.keys():
                if fn in self.pulse_registry.keys():
                    for num in item.pulse_registry[fn]:
                        if num in self.pulse_registry[fn]:
                            return True
        else:
            if item.origfn in self.pulse_registry.keys() and \
                item.number in self.pulse_registry[item.origfn]:
                return True
        return False

    def write_to_file(self, basefn=None):
        """Dump the pulse to file.
            'basefn' is base filename to use. Default will
            be to use the first profile's original data file's 
            name (with the extension).
        """
        #
        # Assume original filename has an extension
        # NOTE: This is probably not good to assume.
        #
        if basefn is None:
            basefn, extension = os.path.splitext(self.origfn)
        file = open("%s.summedprof" % basefn, 'w')
        file.write("# Original data file              = %s\n" % self.origfn)
        file.write("# Pulse Number                    = %d\n" % self.number)
        file.write("# MJD of start of pulse           = %0.15f\n" % self.mjd)
        file.write("# Time into observation (seconds) = %f\n" % self.time)
        file.write("# Duration of pulse (seconds)     = %0.15f\n" % self.duration)
        file.write("# Profile bins                    = %d\n" % self.N)
        file.write("# Width of profile bin (seconds)  = %g\n" % self.dt)
        if self.on_pulse is not None:
            for i, (lo,hi) in enumerate(self.on_pulse):
                file.write("# On-pulse region %2d (phase)      = %f-%f\n" % \
                                                                    (i,lo,hi))
        file.write("# Number of profiles summed       = %d\n" % self.count)
        for fn in self.pulse_registry.keys():
            for num in sorted(self.pulse_registry[fn]):
                file.write("# Pulse registry                  = %s:%d\n" % \
                                                                    (fn, num))
        file.write("###################################\n")
        for i, val in enumerate(self.profile):
            file.write("%-10d %f\n" % (i, val))
        file.close()

def read_pulse_from_file(filename):
    """Read pulse information from 'filename' and
        return a Pulse object.
    """
    file = open(filename, 'r')
    profile = []
    on_pulse_regions = []
    for line in file.readlines():
        if line.startswith("# Original data file"):
            origfn = line.split('=')[-1].strip()
        elif line.startswith("# Pulse Number"):
            number = int(line.split('=')[-1].strip())
        elif line.startswith("# MJD of start of pulse"):
            mjd = float(line.split('=')[-1].strip())
        elif line.startswith("# Time into observation (seconds)"):
            time = float(line.split('=')[-1].strip())
        elif line.startswith("# Duration of pulse (seconds)"):
            duration = float(line.split('=')[-1].strip())
        elif line.startswith("# Width of profile bin (seconds)"):
            dt = float(line.split('=')[-1].strip())
        elif line.startswith("# On-pulse region"):
            lo = float(line.split('=')[-1].split('-')[0].strip())
            hi = float(line.split('=')[-1].split('-')[1].strip())
            on_pulse_regions.append((lo,hi))
        elif line.startswith("#"):
            pass
        else:
            # Profile value
            profile.append(float(line.split()[-1].strip()))
    return Pulse(number, mjd, time, duration, np.array(profile), \
                    origfn, dt, on_pulse_regions)


class OnPulseRegionError(Exception):
    """Error when on-pulse region is ill-defined.
    """
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return "On-pulse region is ill-defined. %s" % self.message
