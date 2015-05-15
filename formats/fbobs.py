"""
A module for reading filterbank observations 
made up of multiple filterbank files.

NOTE: Uses PRESTO's sigproc.py module to read the
      filterbank files' headers.
      Also uses filterbank.py objects to keep track
      of and read the files.

Patrick Lazarus, May 9th, 2010
"""

import numpy as np
import filterbank

class fbobs:
    """A filterbank observation class."""
    def __init__(self, filfns):
        fbs = [filterbank.filterbank(fn) for fn in filfns]
        startmjds = np.array([fb.header['tstart'] for fb in fbs])
        sortinds = startmjds.argsort()

        self.fbs = [fbs[ii] for ii in sortinds]
        self.filenames = [fb.filename for fb in fbs]
        self.numfiles = len(self.fbs)
        self.startmjds = startmjds[sortinds]
        self.nsamps = np.array([fb.number_of_samples for fb in self.fbs])
        self.tsamp = fbs[0].header['tsamp'] # Assume sample time is uniform over observation
        self.lengths = self.nsamps*self.tsamp # length of each file in seconds

        self.endsamps = np.cumsum(self.nsamps)
        self.endtimes = self.endsamps*self.tsamp
        self.startsamps = np.concatenate(([0],self.endsamps[:-1]))
        self.starttimes = self.startsamps*self.tsamp

        self.obslen = self.endtimes[-1]
        self.number_of_samples = self.endsamps[-1]

        self.fbs[0].compute_frequencies()
        self.frequencies = self.fbs[0].frequencies
        self.nchans = self.fbs[0].header['nchans']

    def close_all(self):
        """Close all filterbank files.
        """
        for fb in self.fbs:
            fb.close()

    def __del__(self):
        self.close_all()

    def get_time_interval(self, starttime, endtime):
        """Get a chunk of filterbank observation starting
            at time 'starttime' and ending at time 'endtime'.
            Both arguments are in seconds.
        """
        # Convert to samples and use self.get_sample_interval
        startsamp = starttime/self.tsamp
        endtsamp = endtime/self.tsamp
        return self.get_sample_interval(startsamp, endsamp)

    def get_sample_interval(self, startsamp, endsamp):
        """Get a chunk of filterbank observation starting
            at sample 'startsamp' and ending at sample 'endtime'.
        """
        if startsamp > endsamp:
            raise ValueError("Start of interval must preceed end of interval!")
        if startsamp <= 0:
            startsamp = 0
        if endsamp >= self.endsamps[-1]:
            endsamp = self.endsamps[-1] 
        # Find files at start and end of interval
        for ii, (start, end) in enumerate(zip(self.startsamps, self.endsamps)):
            if start <= startsamp <= end:
                startfile = ii
            if start <= endsamp <= end:
                endfile = ii
        # Need to be careful about intervals spanning multiple files.
        chunks = []
        if startfile == endfile:
            # Grab interval
            firstsamp = startsamp - self.startsamps[startfile]
            numsamps = endsamp - startsamp
            self.fbs[startfile].seek_to_sample(firstsamp)
            chunks.append(self.fbs[startfile].read_Nsamples(numsamps))
        else:
            # Grab start chunk, end chunk and full files in between
            # Be sure to seek to start of data
            firstsamp = startsamp - self.startsamps[startfile]
            numsamps = self.endsamps[startfile] - startsamp
            self.fbs[startfile].seek_to_sample(firstsamp)
            chunks.append(self.fbs[startfile].read_Nsamples(numsamps))
            
            for ii in range(startfile+1,endfile):
                chunks.append(self.fbs[ii].read_all_samples())

            numsamps = endsamp - self.startsamps[endfile]
            self.fbs[endfile].seek_to_data_start()
            chunks.append(self.fbs[endfile].read_Nsamples(numsamps))

        return np.concatenate(chunks)
