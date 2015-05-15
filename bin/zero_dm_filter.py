#!/usr/bin/env python

"""
zero_dm_filter.py - Perform Zero-DM Filter on data.

Patrick Lazarus, Jan 9, 2010
"""
import sys
import optparse
import struct
import warnings

import numpy as np

import sigproc

import filterbank
import colour


def write_header(outfbfile, header_params, header):
    """Given list of header flags and values and an output
        file object write the new filterbank header.
    """
    for hp in header_params:
        hdrval = sigproc.addto_hdr(hp, header[hp])
        outfbfile.write(hdrval)


def filter(data):
    """Perform Zero-DM Filter on data.
        Subtract mean of each time sample from itself.
    """
    
    data_avg = data.mean()
    if data_avg.dtype != data.dtype:
        # Round and cast data_avg to retain data's dtype
        data_avg = np.round(data_avg).astype(data.dtype)
    return data - data_avg


def write_data(infbfile, outfbfile):
    # Loop over file.
    for i in range(0, infbfile.number_of_samples):
        # Read data
        data = infbfile.read_sample()
        # Filter Data
        data = filter(data)
        # Write Data
        data.tofile(outfbfile)


def main():
    sys.stdout.write("Working...")
    sys.stdout.flush()

    # Open infile
    infbfile = filterbank.filterbank(options.infile)

    # Open outfile
    outfbfile = open(options.outname, 'wb')

    # Write header
    header = infbfile.header
    header_params = infbfile.header_params
    write_header(outfbfile, header_params, header)
    sys.stdout.write("\rHeader written... Writing data...")
    sys.stdout.flush()

    # Write file
    write_data(infbfile, outfbfile)
    sys.stdout.write("\rDone!" + " "*50 + "\n")
    sys.stdout.flush()
    
    # Close files
    outfbfile.close()
    infbfile.close()


if __name__=='__main__':
    parser = optparse.OptionParser(usage='%prog [options] infile', \
                description="Perfom Zero-DM Filter on filterbank file.")
    parser.add_option('-o', '--outname', dest='outname', type='string', \
                help="Output filename.", default=None)
    parser.add_option('-d', '--debug', dest='debug', action='store_true', \
                help="Print useful debugging information. " \
                     "(Default: Don't print debug info.)", default=False)
    (options, args) = parser.parse_args()

    if len(args)==0:
        parser.print_help()
        sys.exit(1)
    elif len(args)!=1:
        sys.stderr.write("Only one input file must be provided!\n")
    else:
        options.infile = args[-1]

    if options.outname is None:
        sys.stderr.write("An outname must be provided. " \
                            "(Use -o/--outname on command line).\n")
        sys.exit(1)
    main()
    
