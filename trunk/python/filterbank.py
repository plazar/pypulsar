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
            print "ERROR: File does not exist!\n\t(%s)" % filfn
            self = None
            return
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
            if self.header['nbits'] == 32:
                self.dtype = 'float32'
            else:
                self.dtype = 'int%d' % self.header['nbits']
            self.header_size = self.filfile.tell()
            self.data_size = os.stat(self.filename)[6] - self.header_size
            bytes_per_sample = self.header['nchans'] * (self.header['nbits']/8)
            if self.data_size % bytes_per_sample:
                warnings.warn("Not an integer number of samples in file.")
            self.number_of_samples = self.data_size / bytes_per_sample

    def read_sample(self):
        """
        Read one sample starting at current
        position. Return as numpy array.
        
        NOTE: No checks are made to see if 
              current position is start of
              a sample.
        """
        self.read_header()
        return np.fromfile(self.filfile, dtype=self.dtype, count=self.header['nchans'])

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
        return np.fromfile(self.filfile, dtype=self.dtype, count=self.header['nchans']*N)
    
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
        self.filfile.seek(self.header_size + self.header['nbits']/8*self.header['nchans']*sampnum)
        
    def seek_to_position(self, posn):
        """
        See to position 'posn' relative to
        beginning of file.
        """
        self.filfile.seek(posn)
