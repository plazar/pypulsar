#!/usr/bin/env python

"""
combinefil.py - Combine filterbank files.

Patrick Lazarus, Nov. 27, 2009
"""
import sys
import optparse
import struct
import warnings

import numpy as np

import sigproc

import filterbank
import colour

SAMPLES_PER_READ = 256


def sort_fb_files(fbfiles):
    """Given a list of filterbank objects sort by frequency.
    """
    freqbands = []
    reversed = []
    for fb in fbfiles:
        chanbw = fb.header['foff']
        if chanbw < 0:
            reversed.append(True)
        else:
            reversed.append(False)
        totalbw = chanbw*fb.header['nchans']
        lofreq = fb.header['fch1']
        hifreq = lofreq + totalbw
        freqbands.append([lofreq, hifreq])
    
    reversed = np.array(reversed)
    if not (reversed.all() or np.bitwise_not(reversed)):
        raise ValueError(colour.cstring("Frequency bands are not ordered" \
                                        "the same.", 'error'))
        
    freqbands = np.array(freqbands).astype(float)
    order = freqbands.argsort(axis=0).transpose()[0]
    if reversed.all():
        # reverse order
        if options.debug:
            colour.cprint("Frequency bands are all inverted.", 'debug')
        order = order[::-1]
    freqbands = freqbands[order]
    if np.any(freqbands.flat != sorted(freqbands.flat, reverse=reversed.all())):
        if options.debug:
            colour.cprint(freqbands, 'debug')
        raise ValueError(colour.cstring("Frequency bands have overlaps or " \
                                        "are inverted.", 'error'))
    
    sorted_fbfiles = []
    for i in order:
        sorted_fbfiles.append(fbfiles[i])
    return sorted_fbfiles
  

def write_header(infbfiles, outfbfile):
    """Given list of input filterbank objects and an output
        file object write the new filterbank header.
    """
    infbfiles = sort_fb_files(infbfiles)
    for hp in infbfiles[0].header_params:
        if hp == 'nchans':
            val = sum([f.header['nchans'] for f in infbfiles])
        else:
            val = infbfiles[0].header[hp]
        hdrval = sigproc.addto_hdr(hp, val)
        outfbfile.write(hdrval)
  

def write_data(infbfiles, outfbfile):
    """Given list of input filterbank objects and an output
        file object write the new filterbank data.
    """
    infbfiles = sort_fb_files(infbfiles)
    for i in range(0, infbfiles[0].number_of_samples/SAMPLES_PER_READ):
        data = [fb.read_Nsamples(SAMPLES_PER_READ) for fb in infbfiles]
        for d,fb in zip(data,infbfiles):
            d.shape = (SAMPLES_PER_READ, fb.header['nchans'])
        data = np.atleast_2d(np.hstack(data))
        data.tofile(outfbfile)
    
    # Write leftover samples to outfiles
    num_samples_leftover = infbfiles[0].number_of_samples % SAMPLES_PER_READ
    if num_samples_leftover:
        data = [fb.read_Nsamples(num_samples_leftover) for fb in infbfiles]
        for d,fb in zip(data,infbfiles):
            d.shape = (num_samples_leftover, fb.header['nchans'])
        data = np.atleast_2d(np.hstack(data))
        data.tofile(outfbfile)
        

def main():
    warnings.warn(colour.cstring("Not checking if .fil files are contiguous " \
                                "frequency bands, or the same length, etc.", \
                                'warning'))
    
    sys.stdout.write("Working...")
    sys.stdout.flush()
    
    # Open infiles
    infbfiles = [filterbank.filterbank(infile) for infile in options.infiles]

    outfile = open(options.outname, 'wb')
    
    # Write header
    write_header(infbfiles, outfile)
    sys.stdout.write("\rHeader written... Writing data...")
    sys.stdout.flush()
    
    # Write data
    write_data(infbfiles, outfile)
    sys.stdout.write("\rDone!" + " "*50 + "\n")
    sys.stdout.flush()
    
    outfile.close()
    
    # Close infiles
    for fb in infbfiles:
        fb.close()


if __name__=='__main__':
    parser = optparse.OptionParser(usage='%prog [options] infiles', \
                description="Combine filterbank data files for contiguous " \
                            "frequency bands into a single file.")
    parser.add_option('-o', '--outname', dest='outname', type='string', \
                help="Output filename.", default=None)
    parser.add_option('-d', '--debug', dest='debug', action='store_true', \
                help="Print useful debugging information. " \
                     "(Default: Don't print debug info.)", default=False)
    (options, args) = parser.parse_args()

    if len(args)==0:
        parser.print_help()
        sys.exit(1)
    else:
        options.infiles = args

    if options.outname is None:
        sys.stderr.write("An outname must be provided. " \
                            "(Use -o/--outname on command line).\n")
        sys.exit(1)
    main()
