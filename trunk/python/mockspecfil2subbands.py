#!/usr/bin/env python
#
# Simple script to convert sigproc filterbank data
# into subbands readable by PRESTO. 
#
#               Patrick Lazarus, June 11, 2009

import sys
import os.path
import struct
import optparse
import numpy as np

import sigproc

import filterbank
import coordconv

# Constants
SAMPLES_PER_READ = 1024*4 # Number of samples to read at once

parser = optparse.OptionParser(usage='%prog [options] infile', description="Convert filterbank data (from MockSpec data) to PRESTO subbands. Each subband is one channel.")
parser.add_option('-o', '--outname', dest='outname', type='string', help="Output filename. Do not include extension.", default=None)
(options, args) = parser.parse_args()

if len(args)==0:
    parser.print_help()
    sys.exit(1)
else:
    options.infile = args[0]

if options.outname is None:
    sys.stderr.write("An outname must be provided. (Use -o/--outname on commandline).\n")
    sys.exit(1)

def writeInfoFile(filfile):
    """
    Write .inf files for subbands generated 
    from 'filfile'.
    """
    basefilename = "%s.sub" % options.outname
    inffilename = "%s.inf" % basefilename
    inffile = open(inffilename, 'w')
    inffile.write(" Data file name without suffix          =  %s\n" % \
                    basefilename)
                    
    telescope_id = filfile.header['telescope_id']
    if telescope_id in sigproc.ids_to_telescope.keys():
        telescope = sigproc.ids_to_telescope[telescope_id]
    else:
        telescope = "????"
    inffile.write(" Telescope used                         =  %s\n" % \
                    telescope)
                    
    machine_id = filfile.header['machine_id']
    if machine_id in sigproc.ids_to_machine.keys():
        machine = sigproc.ids_to_machine[machine_id]
    else:
        machine = "????"
    inffile.write(" Instrument used                        =  %s\n" % \
                    machine)

    inffile.write(" Object being observed                  =  %s\n" % \
                    filfile.header['source_name'])

    inffile.write(" J2000 Right Ascension (hh:mm:ss.ssss)  =  %s\n" % \
                    coordconv.rastr_to_fmrastr(filfile.header['src_raj']))

    inffile.write(" J2000 Declination     (dd:mm:ss.ssss)  =  %s\n" % \
                    coordconv.decstr_to_fmdecstr(filfile.header['src_dej']))

    inffile.write(" Data observed by                       =  Unknown\n")

    inffile.write(" Epoch of observation (MJD)             =  %.15f\n" % \
                    filfile.header['tstart'])

    inffile.write(" Barycentered?           (1=yes, 0=no)  =  0\n")

    inffile.write(" Number of bins in the time series      =  %d\n" % \
                    filfile.number_of_samples)

    inffile.write(" Width of each time series bin (sec)    =  %g\n" % \
                    filfile.header['tsamp'])

    inffile.write(" Any breaks in the data? (1=yes, 0=no)  =  0\n")

    inffile.write(" Type of observation (EM band)          =  Radio\n")

    inffile.write(" Beam diameter (arcsec)                 =  175\n") # ALFA

    inffile.write(" Dispersion measure (cm-3 pc)           =  0\n")

    if filfile.header['foff'] < 0:
        chanbw = - filfile.header['foff']
        totalbw = chanbw*filfile.header['nchans']
        lofreq = filfile.header['fch1']-totalbw
    else:
        chanbw = filfile.header['foff']
        totalbw = chanbw*filfile.header['nchans']
        lofreq = filfile.header['fch1']
        
    inffile.write(" Central freq of low channel (Mhz)      =  %f\n" % \
                    lofreq)

    inffile.write(" Total bandwidth (Mhz)                  =  %f\n" % \
                    totalbw)

    inffile.write(" Number of channels                     =  %d\n" % \
                    filfile.header['nchans'])

    inffile.write(" Channel bandwidth (Mhz)                =  %f\n" % \
                    chanbw)

    inffile.write(" Data analyzed by                       =  Patrick Lazarus\n")

    inffile.write(" Any additional notes:\n")
    inffile.write("    Input filterbank file created from MockSpec data (AO)\n")
    inffile.write("    using psrfits2fil, written by Julia Deneva (?)\n")
    inffile.write("    Subbands and inf file created by mockspecfil2subbands.py\n")
    inffile.write("    written by Patrick Lazarus, June 11, 2009\n")

    inffile.close()

def main():
    # Open infile
    infbfile = filterbank.filterbank(options.infile)

    # Write .inf file (open/write/close)
    writeInfoFile(infbfile)

    # Open outfiles
    if infbfile.header['foff'] > 0:
        filenames = ["%s.sub%04d" % (options.outname, subnum) for subnum in range(0,infbfile.header['nchans'])]
    elif infbfile.header['foff'] < 0:
        # Invert subband files with respect to filterbank file.
        filenames = ["%s.sub%04d" % (options.outname, subnum) for subnum in range(infbfile.header['nchans']-1,0-1,-1)]
    else:
        sys.stderr.write("Channel bandwidth is 0! Exiting...\n")
        sys.exit(1)
        
    outfiles = [open(fn, 'wb') for fn in filenames]

    sys.stdout.write("Working...")
    # Loop
    for i in range(1, infbfile.number_of_samples/SAMPLES_PER_READ):
        # Read samples
        data = infbfile.read_Nsamples(SAMPLES_PER_READ)
        # Re-organize data
        data.shape = (SAMPLES_PER_READ, infbfile.header['nchans'])
        data = data.transpose()
        # Loop
        for j in range(0, infbfile.header['nchans']):
            # Write to outfiles
            data[j].tofile(outfiles[j])
    
    # Write leftover samples to outfiles
    num_samples_leftover = infbfile.number_of_samples % SAMPLES_PER_READ
    data = infbfile.read_Nsamples(num_samples_leftover)
    # Re-organize data
    data.shape = (num_samples_leftover, infbfile.header['nchans'])
    data = data.transpose()
    # Loop
    for j in range(0, infbfile.header['nchans']):
        # Write to outfiles
        data[j].tofile(outfiles[j])
    
    # Close infile
    infbfile.close()
    # Close outfiles
    for outfile in outfiles:
        outfile.close()
    
    sys.stdout.write('\rDone!       \n')
    sys.stdout.flush()

if __name__=='__main__':
    main()
