#!/usr/bin/env python

"""
sum_profs.py

Sum profiles.
Command line arguments are used to:
    - locate profiles
    - define post-sum processing on the final profile.

Patrick Lazarus, June 30, 2009
"""

import sys
import glob
import optparse
import types
import pulse

post_sum_processing = []

def main():
    pulsefiles = args + glob.glob(options.glob_expr)
    if len(pulsefiles) < 2:
        print pulsefiles
        sys.stderr.write("Only %d pulse files provided. Exiting!\n" % \
                                len(pulsefiles))
        sys.exit(1)

    print "Summing %d profiles" % len(pulsefiles)
   
    psum = pulse.read_pulse_from_file(pulsefiles[0]) + \
            pulse.read_pulse_from_file(pulsefiles[1])
    for pulsefile in pulsefiles[2:]:
        psum += pulse.read_pulse_from_file(pulsefile)

    # Do post processing
    for procstep in post_sum_processing:
        method_name, method_args = procstep
        method = eval('psum.%s' % method_name)
        if type(method_args)==types.TupleType:
            print "Applying %s with arguments %s" % (method_name, method_args)
            method(*method_args)
        elif method_args is None:
            print "Applying %s" % method_name
            method()
        else:
            print "Applying %s with argument %s" % (method_name, method_args)
            method(method_args)

    # Write to file
    psum.write_to_file(basefn=options.outname)


def register_post_sum_processing(option, opt_str, value, parser):
    """Register a post-sum processing step.
    """
    post_sum_processing.append((option.dest, value))


if __name__=='__main__':
    parser = optparse.OptionParser(usage="%prog [options] infiles", description="Given Pulse object files sum them. Optionally apply processing to the summed profile. Output summed profile to file. Written by Patrick Lazarus.", version="%prog v0.9 (by Patrick Lazarus)", prog="sum_profs.py")
    parser.add_option('-g', '--glob-expr', dest='glob_expr', type='string', help="Glob expression identifying prof files to sum. (Default: Nothing).", default="")
    parser.add_option('-o', '--outname', dest='outname', type='string', help="Base filename of output summed-profile. (Default: use original datafile's basename. If pulses come from multiple data file then one of them is picked).", default=None) 

    postsum_group = optparse.OptionGroup(parser, "Post-Sum Processing", "The following options define post-sum processing to be applied to the summed profile before writing it to file. The order in which the options are provided on the command line is important.")
    postsum_group.add_option('--scale', dest='scale', nargs=0, action='callback', callback=register_post_sum_processing, help="Scale the profile by subtracting the mean of the off-pulse region and dividing by the stddev of the off-pulse region.")
    postsum_group.add_option('--downsample', dest='downsample', type='int', action='callback', callback=register_post_sum_processing, help="Downsample profile by given factor.")
    postsum_group.add_option('--smooth', dest='smooth', type='int', action='callback', callback=register_post_sum_processing, help="Smooth profile using a boxcar of given width.")
    postsum_group.add_option('--detrend', dest='detrend', type='int', action='callback', callback=register_post_sum_processing, help="Break profile into N chunks (where N is provided on command line) and subtract a linear trend from each chunk.")
    postsum_group.add_option('--interpolate', dest='interpolate', type='int', action='callback', callback=register_post_sum_processing, help="Interpolate profile so it has N bins across it (where N is provided on the command line).")
    parser.add_option_group(postsum_group)
    
    options, args = parser.parse_args()
    main()
