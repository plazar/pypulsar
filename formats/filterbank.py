#!/usr/bin/env python
#
# A module for reading filterbank files.
# 
# NOTE: Uses PRESTO's sigproc.py module to read the
#       filterbank file's header
#
#           Patrick Lazarus, June 6th, 2009

import sys
import warnings
import os.path
import numpy as np
import sigproc

class filterbank:
    def __init__(self, filfn):
        if not os.path.isfile(filfn):
            raise ValueError("ERROR: File does not exist!\n\t(%s)" % filfn)
        else:
            self.filename = filfn
            self.already_read_header = False
            self.header_params = []
            self.header = {}
            self.header_size = None # size of header in bytes
            self.data_size = None
            self.number_of_samples = None
            self.dtype = None
            self.filfile = open(filfn, 'rb')
            self.read_header()
            self.compute_frequencies()

    def close(self):
        self.filfile.close()

    def __del__(self):
        self.close()

    def read_header(self):
        if not self.already_read_header:
            self.header_params = []
            self.header = {}
            self.already_read_header = True
            self.seek_to_header_start()
            param = ""
            while (param != 'HEADER_END'):
                param, val = sigproc.read_hdr_val(self.filfile)
                self.header[param] = val
                self.header_params.append(param)
            
            # Calculate additional information
            # Such as: datatype, numsamps, datasize, hdrsize
            if self.nbits == 32:
                self.dtype = 'float32'
            else:
                self.dtype = 'uint%d' % self.nbits
            self.header_size = self.filfile.tell()
            self.data_size = os.stat(self.filename)[6] - self.header_size
            bytes_per_sample = self.nchans * (self.nbits/8)
            if self.data_size % bytes_per_sample:
                warnings.warn("Not an integer number of samples in file.")
            self.number_of_samples = self.data_size / bytes_per_sample

    def __getattr__(self, name):
        return self.header[name]

    def print_header(self):
        """Print header parameters and values.
        """
        self.read_header()
        for param in self.header_params:
            print "%s: %s" % (param, self.header[param])

    def compute_frequencies(self):
        """Compute frequencies (in MHz) and store in 
            attribute self.frequencies.

            NOTE: frequencies correspond to order of channels
                    in filterbank file.
        """
        self.read_header()
        self.frequencies = self.fch1 + self.foff*np.arange(self.nchans)
        self.is_hifreq_first = (self.foff < 0)

    def read_sample(self):
        """
        Read one sample starting at current
        position. Return as numpy array.
        
        NOTE: No checks are made to see if 
              current position is start of
              a sample.
        """
        self.read_header()
        return np.fromfile(self.filfile, dtype=self.dtype, count=self.nchans)

    def read_all_samples(self):
        """
        Read all samples from file.
        Return as numpy array.
        """
        self.seek_to_data_start()
        return np.fromfile(self.filfile, dtype=self.dtype)

    def read_Nsamples(self, N):
        """
        Read N samples starting at current
        position. Return as numpy array.
        
        NOTE: No checks are made to see if 
              current position is start of
              a sample.
        """
        self.read_header()
        return np.fromfile(self.filfile, dtype=self.dtype, count=self.nchans*N)
    
    def seek_to_header_start(self):
        self.filfile.seek(0)

    def seek_to_data_start(self):
        self.read_header()
        self.filfile.seek(self.header_size)

    def seek_to_sample(self, sampnum):
        """
        Seek to sample 'sampnum'. First sample is
        sampnum=0.
        """
        self.read_header()
        self.filfile.seek(self.header_size + self.nbits/8*self.nchans*sampnum)
        
    def seek_to_position(self, posn):
        """
        See to position 'posn' relative to
        beginning of file.
        """
        self.filfile.seek(posn)
