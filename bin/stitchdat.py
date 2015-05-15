#!/usr/bin/env python
"""
Stitch together .dat file to form a longer observation.

Patrick Lazarus, May 30th, 2010
"""
import optparse
import os.path
import sys
import warnings
import copy
import numpy as np
import psr_utils
from pypulsar.formats import datfile


def sort_by_mjd(datfiles):
    """Given a list of Datfile objects sort it in place by start time.
    """
    cmp_mjd = lambda dat1, dat2: cmp(dat1.infdata.epoch, dat2.infdata.epoch)
    datfiles.sort(cmp=cmp_mjd)


def main():
    warnings.warn("Not checking if all .dat files have same observing band and sample time.")
    datfiles = [datfile.Datfile(datfn) for datfn in options.infiles]
    sort_by_mjd(datfiles)

    # Open output file
    outfile = open(options.outname+".dat", 'wb')
    numsamps = 0
    print "Working on", os.path.split(datfiles[0].datfn)[1]
    dtype = datfiles[0].dtype
    data = datfiles[0].read_all()
    datfiles[0].close()
    data.tofile(outfile)
    numsamps += data.size
    for ii in range(1,len(datfiles)):
        print "Working on", os.path.split(datfiles[ii].datfn)[1]
        # Calculate number of padvalues required
        mjd_diff = datfiles[ii].currmjd_actual - datfiles[ii-1].currmjd_actual
        sec_diff = mjd_diff * psr_utils.SECPERDAY
        samp_diff = sec_diff / datfiles[ii].infdata.dt
        warnings.warn("Padding by integer number of sample times.")
        numpadvals = int(np.around(samp_diff))

        # Generate array for padding
        padval = np.median(data)
        padvals = np.ones(numpadvals, dtype=dtype)*padval
        if options.debug:
            print "Padding by %d samples" % numpadvals
            print "Value used for padding: %g" % padval
            print "Padding by integer number of bins caused %f bins to be " \
                    "discarded/added" % (samp_diff-numpadvals)
        padvals.tofile(outfile)
        numsamps += padvals.size

        # Write next dataset
        data = datfiles[ii].read_all()
        datfiles[ii].close()
        data.tofile(outfile)
        numsamps += data.size

    # Close output file
    outfile.close()

    # Generate new .inf file
    inf = copy.deepcopy(datfiles[0].infdata)
    inf.N = numsamps
    inf.tofile(options.outname+'.inf')

    print "Total number of samples written:", numsamps

if __name__=='__main__':
    parser = optparse.OptionParser(usage='%prog [options] infiles', \
                description="Stitch together multiple .dat files to" \
                            "form a longer observation. Padding will" \
                            "be performed as needed.", \
                prog="stitchdat.py", version="Patrick Lazarus, v0.1")
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

    if len(options.infiles) < 2:
        sys.stderr.write("Need at least 2 files to stitch together.\n")
        sys.exit(2)

    if options.outname is None:
        sys.stderr.write("An outname must be provided. " \
                            "(Use -o/--outname on command line).\n")
        sys.exit(3)
    main()
