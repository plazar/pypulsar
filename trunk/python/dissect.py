#!/usr/bin/env python

"""
dissect.py

Cut a timeseries into individual pulses.
* Can be used with:
    - constant period
    - polycos (not implemented yet)
    - parfile (not implemented yet)
* Can dump pulses with SNR above given threshold to:
    - .bestprof files (not implemented yet)
    - .dat files (not implemented yet)

Patrick Lazarus, June 16, 2009
"""

import sys
import optparse
import warnings
import numpy as np
import datfile
import mypolycos

import ppgplot
# Import matplotlib/pylab and set for non-interactive plots
import matplotlib
matplotlib.use('Agg')
import pylab as plt

# Constants
JOYDIV_SEP = 0.5    # vertical separation per profile in same units 
                    # profiles are plotted in for joydiv plot.
DEFAULT_WIDTHS = [1,2,4,8,16,32] # Powers of 2
# DEFAULT_WIDTHS = [1,2,3,4,6,9,14,20,30] # Default from single_pulse_search.py
                    # default boxcar widths to use for smoothing
                    # when searching for pulses.

   

def main():
    # Open file
    datfn = args[0]
    timeseries = datfile.Datfile(datfn)

    if options.parfile is not None:
        # generate polycos
        # get periods from polycos
        polycos = mypolycos.create_polycos(options.parfile, timeseries.infdata)
        get_period = lambda mjd: 1.0/polycos.get_phs_and_freq(int(mjd), \
                                                               mjd-int(mjd))[1]
    elif options.polycofile is not None:
        raise NotImplementedError("--use-polycos option in dissect.py is not implemented yet")
    elif options.period is not None:
        get_period = lambda mjd: options.period
    else:
        raise "Unknown option for reading periods!"

    # Loop over pulses in timeseries. Examine pulses one at a time.
    good_pulses = []
    snrs = []
    notes = []
    for current_pulse in timeseries.pulses(get_period):
        maxsnr = 0
        bestfactor = 0
        current_pulse.set_onoff_pulse_regions([(options.on_pulse_start, \
                                                options.on_pulse_end)])
        if current_pulse.is_masked(numchunks=5): 
            continue
        for numbins in options.widths:
            pulse = current_pulse.make_copy()
            pulse.smooth(numbins)
            snr = get_snr(pulse)
            if snr > options.threshold:
                if snr > maxsnr:
                    if maxsnr==0 and bestfactor==0:
                        # First time snr is above threshold
                        snrs.append(snr)
                        notes.append("smoothed by %3d bins" % numbins)
                        good_pulses.append(current_pulse)
                    else:
                        # Better smootfactor/snr found, update, but don't
                        # add pulse to good_pulses again.
                        snrs[-1] = snr
                        notes[-1] = "smoothed by %3d bins" % numbins
                    maxsnr = snr
                    bestfactor = numbins
            

    print_report(good_pulses, snrs=snrs, notes=notes)
    if options.create_output_files:
        if options.create_text_files:
            write_pulses(good_pulses, timeseries)
        if options.create_plot_files:
            plot_pulses(good_pulses, timeseries, options.downfactor)
        if options.create_joydiv_plot:
            #joy_division_plot2(good_pulses, timeseries.pulses(get_period), \
            #                timeseries, options.downfactor, \
            #                options.on_pulse_start, options.on_pulse_end)
            joy_division_plot(good_pulses, timeseries, options.downfactor)


def get_snr(pulse, uncertainty=1):
    """Compute and return statistics about the given pulse.
        Also accepts uncertainty on profile bins. 'uncertainty'
        can be a scalar, or numpy array.
    """
    # Scale profile
    copy_of_pulse = pulse.make_copy()
    copy_of_pulse.scale()
    snr = np.max(copy_of_pulse.get_on_pulse())
    # snr = np.sum(copy_of_pulse.get_on_pulse())
    warnings.warn("Only checking on-pulse region for pulses.")
    return snr


def print_report(pulses, snrs=None, notes=None):
    """Print a report given the pulses provided.
    """
    print "Autopsy report:"
    print "\tNumber of good pulses found: %s" % len(pulses)
    use_snrs = ""
    use_notes = ""
    if snrs is not None and len(snrs) == len(pulses):
        use_snrs = "SNR"
    if notes is not None and len(notes) == len(pulses):
        use_notes = "Notes"
    print "%s%s%s%s%s%s" % ("#".center(7), "MJD".center(15), \
                                "Time".center(11), "Duration".center(13), \
                                use_snrs.center(9), use_notes)
    for i, pulse in enumerate(pulses):
        sys.stdout.write(("%d" % pulse.number).center(7))
        sys.stdout.write(("%5.4f" % pulse.mjd).center(15))
        sys.stdout.write(("%5.2f" % pulse.time).center(11))
        sys.stdout.write(("%2.4f" % pulse.duration).center(13))
        if use_snrs:
            sys.stdout.write(("%4.2f" % snrs[i]).center(9))
        if use_notes:
            sys.stdout.write("%s" % notes[i])
        sys.stdout.write("\n")
   

def plot_pulses(pulses, timeseries, downfactor=1):
    """Plot each pulse into a separate file.
        Downsample profiles by factor 'downfactor' before plotting.
    """
    for pulse in pulses:
        copy_of_pulse = pulse.make_copy()
        copy_of_pulse.downsample(downfactor)
        plt.figure()
        plt.plot(copy_of_pulse.profile, 'k-', lw=0.5)
        plt.xlabel("Profile bin")
        plt.title("Pulse #%d" % pulse.number)
        plt.savefig("%s.prof%d.ps" % (timeseries.basefn, pulse.number), \
                        orientation='landscape')


def joy_division_plot(pulses, timeseries, downfactor=1):
    """Plot each pulse profile on the same plot separated
        slightly on the vertical axis.
        'timeseries' is the Datfile object dissected.
        Downsample profiles by factor 'downfactor' before plotting.
    """
    first = True
    ppgplot.pgbeg("%s.joydiv.ps/CPS" % timeseries.basefn, 1, 1)
    ppgplot.pgpap(10.25, 8.5/11.0) # Letter landscape
    # ppgplot.pgpap(7.5, 11.7/8.3) # A4 portrait, doesn't print properly
    ppgplot.pgiden()
    ppgplot.pgsci(1)
    
    # Set up main plot
    ppgplot.pgsvp(0.1, 0.9, 0.1, 0.8)
    ppgplot.pglab("Profile bin", "Single pulse profiles", "")
    to_plot = []
    xmin = 0
    xmax = None
    ymin = None
    ymax = None
    for pulse in pulses:
        vertical_offset = (pulse.number-1)*JOYDIV_SEP
        copy_of_pulse = pulse.make_copy()
        copy_of_pulse.downsample(downfactor)
        copy_of_pulse.scale()
        if first:
            summed_prof = copy_of_pulse.profile.copy()
            first = False
        else:
            summed_prof += copy_of_pulse.profile
        prof = copy_of_pulse.profile + vertical_offset
        min = prof.min()
        if ymin is None or min < ymin:
            ymin = min
        max = prof.max()
        if ymax is None or max > ymax:
            ymax = max
        max = prof.size-1
        if xmax is None or max > xmax:
            xmax = max
        to_plot.append(prof)
    yspace = 0.1*ymax
    ppgplot.pgswin(0, xmax, ymin-yspace, ymax+yspace)
    for prof in to_plot:
        ppgplot.pgline(np.arange(0,prof.size), prof)
    ppgplot.pgbox("BNTS", 0, 0, "BC", 0, 0)

    # Set up summed profile plot
    ppgplot.pgsvp(0.1, 0.9, 0.8, 0.9)
    ppgplot.pglab("", "Summed profile", "Pulses from %s" % timeseries.datfn)
    summed_prof = summed_prof - summed_prof.mean()
    ppgplot.pgswin(0, xmax, summed_prof.min(), summed_prof.max())
    ppgplot.pgline(np.arange(0, summed_prof.size), summed_prof)
    ppgplot.pgbox("C", 0, 0, "BC", 0, 0)
    ppgplot.pgclos()


def joy_division_plot2(good_pulses, all_pulses, timeseries, downfactor=1, \
                        on_pulse_start=0.0, on_pulse_end=1.0):
    """Plot each pulse profile on the same plot separated
        slightly on the vertical axis.
        'good_pulses' and 'all_pulses' are generator/list that
            returns/contains Pulse objects. 
        'timeseries' is the Datfile object dissected.
        Downsample profiles by factor 'downfactor' before plotting.
        'on_pulse_start' and 'on_pulse_end' define the on-pulse region
        in terms of rotational phase (between 0.0 and 1.0).

        ** Same as joy_division_plot, but sum all profiles between good
            profiles and display with dotted lines. **
    """
    first = True
    ppgplot.pgbeg("%s.joydiv2.ps/CPS" % timeseries.basefn, 1, 1)
    ppgplot.pgpap(10.25, 8.5/11.0)
    ppgplot.pgiden()
    ppgplot.pgsci(1)
    
    # Set up main plot
    ppgplot.pgsvp(0.1, 0.9, 0.1, 0.8)
    ppgplot.pglab("Profile bin", "Single pulse profiles", "")
    to_plot = []
    sums_to_plot = []
    xmin = 0
    xmax = None
    ymin = None
    ymax = None
    # Make a list of good pulse numbers:
    good_nums = [pulse.number for pulse in good_pulses]
    in_between_pulses_prof = None # Sum of profiles where no pulse is detected
    num_sum = 0
    count = 0
    for pulse in all_pulses:
        pulse.set_onoff_pulse_regions([(on_pulse_start, on_pulse_end)])
        copy_of_pulse = pulse.make_copy()
        copy_of_pulse.downsample(downfactor)
        copy_of_pulse.scale()
        if pulse.number in good_nums:
            if in_between_pulses_prof is not None:
                # First take care of the in between good pulses sum
                vertical_offset = (float(num_sum)/count)*JOYDIV_SEP
                in_between_pulses_prof += vertical_offset
                min = in_between_pulses_prof.min()
                if ymin is None or min < ymin:
                    ymin = min
                max = in_between_pulses_prof.max()
                if ymax is None or max > ymax:
                    ymax = max
                max = in_between_pulses_prof.size-1
                if xmax is None or max > xmax:
                    xmax = max
                sums_to_plot.append(in_between_pulses_prof)
                num_sum = 0
                count = 0
                in_between_pulses_prof = None
            # Now take care of the good pulse
            vertical_offset = (pulse.number-1)*JOYDIV_SEP
            if first:
                summed_prof = copy_of_pulse.profile.copy()
                first = False
            else:
                summed_prof += copy_of_pulse.profile
            prof = copy_of_pulse.profile + vertical_offset
            min = prof.min()
            if ymin is None or min < ymin:
                ymin = min
            max = prof.max()
            if ymax is None or max > ymax:
                ymax = max
            max = prof.size-1
            if xmax is None or max > xmax:
                xmax = max
            to_plot.append(prof)
        else:
            num_sum += pulse.number
            count += 1
            prof = copy_of_pulse.profile.copy()
            if in_between_pulses_prof is None:
                in_between_pulses_prof = prof
            else:
                if in_between_pulses_prof.size > prof.size:
                    in_between_pulses_prof = prof + np.resize(in_between_pulses_prof, prof.shape)
                elif in_between_pulses_prof.size < prof.size:
                    in_between_pulses_prof += np.resize(prof, in_between_pulses_prof.shape)
                else:
                    in_between_pulses_prof += prof
    # Plot last bit of data
    if in_between_pulses_prof is not None:
        vertical_offset = (float(num_sum)/count)*JOYDIV_SEP
        in_between_pulses_prof += vertical_offset
        min = in_between_pulses_prof.min()
        if ymin is None or min < ymin:
            ymin = min
        max = in_between_pulses_prof.max()
        if ymax is None or max > ymax:
            ymax = max
        max = in_between_pulses_prof.size-1
        if xmax is None or max > xmax:
            xmax = max
        sums_to_plot.append(in_between_pulses_prof)
        num_sum = 0
        count = 0
        in_between_pulses_prof = None
    yspace = 0.1*ymax
    ppgplot.pgswin(0, xmax, ymin-yspace, ymax+yspace)
    # Plot partial summed pulses
    ppgplot.pgsci(2)
    for prof in sums_to_plot:
        ppgplot.pgline(np.arange(0,prof.size), prof)
    ppgplot.pgsci(1)
    # Plot good pulses
    for prof in to_plot:
        ppgplot.pgline(np.arange(0,prof.size), prof)
    # Plot boundary of on-pulse region
    ppgplot.pgsls(2)
    ppgplot.pgline(np.array([options.on_pulse_start * prof.size, \
                             options.on_pulse_start * prof.size]), \
                             np.array([ymin,ymax]))
    ppgplot.pgline(np.array([options.on_pulse_end * prof.size, \
                             options.on_pulse_end * prof.size]), \
                             np.array([ymin,ymax]))
    ppgplot.pgsls(1)
    ppgplot.pgbox("BNTS", 0, 0, "BCNTS", 0, 0)

    # Set up summed profile plot
    ppgplot.pgsvp(0.1, 0.9, 0.8, 0.9)
    ppgplot.pglab("", "Summed profile", "Pulses from %s" % timeseries.datfn)
    summed_prof = summed_prof - summed_prof.mean()
    ppgplot.pgswin(0, xmax, summed_prof.min(), summed_prof.max())
    ppgplot.pgline(np.arange(0, summed_prof.size), summed_prof)
    ppgplot.pgbox("C", 0, 0, "BC", 0, 0)
    ppgplot.pgclos()


def write_pulses(pulses, timeseries):
    """Dump the pulses provided to files.
        'timeseries' is the Datfile object dissected.
    """
    for pulse in pulses:
        pulse.write_to_file()

def parse_boxcar_widths(option, opt_str, value, parser):
    """Parse list of boxcar widths from command line.
    """
    widths = [int(w) for w in value.split(',')]
    setattr(parser.values, option.dest, widths)


if __name__ == '__main__':
    parser = optparse.OptionParser(usage="%prog --use-parfile PARFILE | --use-polycos POLYCOFILE | -p --period PERIOD [options] infile.dat", description="Given a input timeseries (a PRESTO .dat file) dissect it into individual pulses and record the pulses that surpass the signification threshold.", version="%prog 0.9")
    parser.add_option('-t', '--threshold', dest='threshold', type='float', action='store', help="Only record pulses more significant than this threshold. (Default: 5).", default=5)
    parser.add_option('-n', '--no-output-files', dest='create_output_files', action='store_false', help="Do not create any output file for each significant pulse detected. (Default: create output files).", default=True)
    parser.add_option('--no-text-files', dest='create_text_files', action='store_false', help="Do not create text file for each significant pulse detected. (Default: create text files).", default=True)
    parser.add_option('-s', '--on-pulse-start', dest='on_pulse_start', type='float', help="Define start of on-pulse region. Value should be a float between 0 and 1. (Default: 0.0).", default=0.0)
    parser.add_option('-e', '--on-pulse-end', dest='on_pulse_end', type='float', help="Define end of on-pulse region. Value should be a float between 0 and 1. (Default: 1.0).", default=1.0)
    parser.add_option('-w', '--widths', dest='widths', type='string', action='callback', callback=parse_boxcar_widths, help="Boxcar widths (in number of samples) to use for smoothing profiles when searching for pulses. widths should be comma-separated _without_ spaces. (Default: Smooth with boxcar widths %s)" % DEFAULT_WIDTHS, default=DEFAULT_WIDTHS)

    period_group = optparse.OptionGroup(parser, "Period Determination", "The following options are different methods for determine the spin period of the pulsar. Exactly one of these options must be provided.")
    period_group.add_option('--use-parfile', dest='parfile', type='string', action='store', help="Determine spin period from polycos generated by tempo using provided parfile.", default=None)
    period_group.add_option('--use-polycos', dest='polycofile', type='string', action='store', help="Determine spin period from polycos in the provided file.", default=None)
    period_group.add_option('-p', '--use-period', dest='period', type='float', action='store', help="Use constant period provided (in seconds).", default=None)
    parser.add_option_group(period_group)

    plot_group = optparse.OptionGroup(parser, "Plotting Options", "The following options affect only the output plots, not the data searching.")
    plot_group.add_option('-d', '--downsample', dest='downfactor', type='int', action='store', help="Down sample profiles by this factor before plotting. (Default: Don't downsample).", default=1)
    parser.add_option('--no-pulse-plots', dest='create_plot_files', action='store_false', help="Do not create plots for each significant pulse detected. (Default: create plots).", default=True)
    parser.add_option('--no-joydiv-plot', dest='create_joydiv_plot', action='store_false', help="Do not create Joy Division plot, where every profile is plotted in a single axes separated slightly in the vertical direction. (Default: create JoyDiv plot).", default=True)

    parser.add_option_group(plot_group)

    options, args = parser.parse_args()

    # Count number of period determination options are provided

    if (options.parfile is not None) + \
        (options.polycofile is not None) + \
        (options.period is not None) != 1:
        parser.print_usage()
        sys.stderr.write("Exactly one (1) period determination option must be provided! Exiting...\n\n")
        sys.exit(1)
    main()
