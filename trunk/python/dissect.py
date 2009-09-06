#!/usr/bin/env python

"""
dissect.py

Cut a timeseries into individual pulses.
* Can be used with:
    - constant period
    - polycos
    - parfile
* Can dump pulses with SNR above given threshold to:
    - .prof files (similar to .bestprof files)

Patrick Lazarus, June 16, 2009
"""

import sys
import os.path
import optparse
import warnings
import numpy as np
import psr_utils
import fftfit
import datfile
import mypolycos
import telescopes
import colour

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

    if options.shift_phase != 0.0:
        # Only shift phase by a fraction of one rotation
        options.shift_phase -= int(options.shift_phase)
        if options.shift_phase < 0.0:
            options.shift_phase += 1.0
    else:
        shift_time = 0.0

    print "Searching %s for single pulses." % timeseries.datfn

    if options.parfile is not None:
        print "Using parfile: %s" % options.parfile
        # generate polycos
        print "Automatically generating polycos..."
        polycos = mypolycos.create_polycos(options.parfile, timeseries.infdata)
        mjd = timeseries.infdata.epoch
        mjdi = int(mjd) # integer part of mjd
        mjdf = mjd-mjdi # fractional part of mjd
        phase, freq = polycos.get_phs_and_freq(mjdi, mjdf)
        if options.debug:
            colour.cprint("MJD at start of file: %r" % mjd, 'debug')
            colour.cprint("Phase at start of file: %f" % phase, 'debug')
        if options.shift_phase != 0.0:
            prof_start_phase = options.shift_phase
            shift_phase = options.shift_phase - phase
            if shift_phase < 0.0:
                shift_phase += 1.0
            shift_time = shift_phase * 1.0/freq
        else:
            prof_start_phase = phase
        # get periods from polycos
        get_period = lambda mjd: 1.0/polycos.get_phs_and_freq(int(mjd), \
                                                               mjd-int(mjd))[1]
    elif options.polycofile is not None:
        print "Using polycos file: %s" % options.polycos
        polycos = mypolycos.polycos(options.polycos)
        mjd = timeseries.infdata.epoch
        mjdi = int(mjd) # integer part of mjd
        mjdf = mjd-mjdi # fractional part of mjd
        phase, freq = polycos.get_phs_and_freq(mjdi, mjdf)
        if options.debug:
            colour.cprint("MJD at start of file: %r" % mjd, 'debug')
            colour.cprint("Phase at start of file: %f" % phase, 'debug')
        if options.shift_phase != 0.0:
            prof_start_phase = options.shift_phase
            shift_phase = options.shift_phase - phase
            if shift_phase < 0.0:
                shift_phase += 1.0
            shift_time = shift_phase * 1.0/freq
        else:
            prof_start_phase = phase
        # get periods from polycos
        get_period = lambda mjd: 1.0/polycos.get_phs_and_freq(int(mjd), \
                                                               mjd-int(mjd))[1]
    elif options.period is not None:
        print "Using constant period: %f" % options.period
        if options.shift_phase != 0.0:
            shift_time = options.shift_phase * options.period
        get_period = lambda mjd: options.period
    else:
        raise ValueError("Unknown option for reading periods.")

    print "On-pulse regions will be set to: %s" % \
            ','.join(['%s:%s' % t for t in options.on_pulse_regions])
    print "Boxcar widths to be used: %s" % \
            ', '.join(['%s' % w for w in options.widths])
    print "Single-pulse SNR threshold: %s" % options.threshold

    # Loop over pulses in timeseries. Examine pulses one at a time.
    good_pulses = []
    snrs = []
    notes = []
    nummasked = 0
    numpulses = 0
    for current_pulse in timeseries.pulses(get_period, \
                                    time_to_skip=shift_time):
        numpulses += 1
        maxsnr = 0
        bestfactor = 0
        current_pulse.set_onoff_pulse_regions(options.on_pulse_regions)
        if current_pulse.is_masked(numchunks=5) and not options.no_toss: 
            nummasked += 1
            continue
        for numbins in options.widths:
            pulse = current_pulse.make_copy()
            pulse.smooth(numbins)
            snr = get_snr(pulse)
            if np.isnan(snr) or snr < 0:
                snr = 0
            if snr > options.threshold:
                if snr >= maxsnr:
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

    print_report(good_pulses, numpulses, nummasked, snrs=snrs, notes=notes, \
                    quiet=options.quiet)
    if options.create_output_files and len(good_pulses) > 0:
        if options.create_text_files:
            print "Writing pulse text files..."
            write_pulses(good_pulses, timeseries)
        if options.create_plot_files:
            print "Creating pulse plots..."
            plot_pulses(good_pulses, timeseries, options.downfactor)
        if options.create_joydiv_plot:
            #joy_division_plot2(good_pulses, timeseries.pulses(get_period), \
            #                timeseries, options.downfactor, \
            #                options.on_pulse_start, options.on_pulse_end)
            print "Making JoyDiv plot..."
            joy_division_plot(good_pulses, timeseries, options.downfactor, \
                                        options.heightstretch)


    if (options.polycofile is not None or options.parfile is not None) and \
                options.write_toas and len(good_pulses) > 0:
        numtoas = 0
        print "Generating TOAs. Please wait..."
        print "TOA threshold:", options.toa_threshold
        print "Min number of pulses for a TOA:", options.min_pulses
        print "Profile template used:", options.template
        # Extract second column from template file
        # First column is index
        template = np.loadtxt(options.template, usecols=(1,))
        # Get TOAs and write them to stdout
        current_pulse = None
        numpulses = 0
        for pulse in good_pulses:
            if current_pulse is None:
                current_pulse = pulse.to_summed_pulse()
                numpulses = 1
            else:
                current_pulse += pulse
                numpulses += 1
            if numpulses < options.min_pulses:
                continue
            if get_snr(current_pulse) > options.toa_threshold:
                # Interpolate and downsample current_pulse so
                # it is same size as template profile
                current_pulse.interp_and_downsamp(template.size)
                # current_pulse.downsample_Nbins(template.size) ## ADDED FOR DEBUGGING
                current_pulse.scale()
                pulseshift, templateshift = write_toa(current_pulse, \
                            polycos, template, timeseries, prof_start_phase, \
                            options.debug)
                numtoas += 1
                if options.write_toa_files:
                    # TOA profiles are already downfactored
                    # do not downfactor more when creating plot
                    plot_toa(numtoas, current_pulse, template, \
                            pulseshift, templateshift)
                    current_pulse.write_to_file("TOA%d" % numtoas)
                current_pulse = None
                numpulses = 0
        print "Number of TOAs: %d" % numtoas
        print "Number of pulses thrown out because 'min pulses' requirement " \
                "or SNR threshold not met: %d" % numpulses


def plot_toa(numtoa, pulse, template=None, pulseshift=0, \
                templateshift=0, basefn=""):
    """Plot the profile used for a TOA.
        - 'numtoa' is the TOA's number within an observation.
        - 'pulse' is the SummedPulse object used to generate 
                the TOA.
        - 'template' is the template used for generating the
                TOA. If it is provided it will be overlayed 
                on the plot. 'template' should be a np.array.
        - 'pulseshift' is an amount of phase to shift pulse by.
        - 'templateshift' is an amount of phase template is shifted by.
        - 'basefn' is the base of the output filename to use.
    """
    if basefn:
        outfn = "%s.TOA%d.ps" % (basefn, numtoa)
    else:
        outfn = "TOA%d.ps" % numtoa
    
    # scale pulse
    copy_of_pulse = pulse.make_copy()
    copy_of_pulse.scale()
    phases = np.linspace(0, 1.0, copy_of_pulse.N)
    tshifted_phases = (phases-templateshift+pulseshift) % (1.0+1e-7)
    tsorted_phase_indices = np.argsort(tshifted_phases)
    # Re-order template
    shifted_template = template[tsorted_phase_indices]
    plt.figure()
    plt.plot(phases, copy_of_pulse.profile, 'k-', lw=0.5)
    if template is not None:
        plt.plot(phases, shifted_template, 'k:', lw=0.5)
    plt.xlabel("Phase (%d profile bins)" % copy_of_pulse.N)
    plt.ylabel("SNR")
    plt.title("TOA #%d" % numtoa)
    plt.savefig(outfn, orientation='landscape')


def write_toa(summed_pulse, polycos, template_profile, \
                        timeseries, start_phase=0.0, debug=False):
    """Given a SummedPulse generate a TOA and write it to stdout. 
        A polycos file is required. 'template_profile' is simply 
        a numpy array. 'timeseries' is a Datfile object.
        'start_phase' is the phase at the start of the profile.
        'debug' is a boolean value that determines if debugging
        info should be displayed.
        Returns shift required to line up template and pulse.
    """
    if template_profile is None:
        raise ValueError("A template profile MUST be provided.")
    # This code is taken from Scott Ransom's PRESTO's get_TOAs.py
    mjdi = int(summed_pulse.mjd) # integer part of MJD
    mjdf = summed_pulse.mjd - mjdi # fractional part of MJD
    (phs, freq) = polycos.get_phs_and_freq(mjdi, mjdf)
    phs -= start_phase
    period = 1.0/freq
    
    # Caclulate offset due to shifting channels to account for DM
    # Hifreq doesn't have a half-channel offset 
    # (why? because get_TOAs.py doesn't. Why...)
    # Why subtract 1 channel to get hifreq?
    hifreq = timeseries.infdata.lofreq + timeseries.infdata.BW - \
                timeseries.infdata.chan_width
    midfreq = timeseries.infdata.lofreq - 0.5*timeseries.infdata.chan_width + \
                0.5*timeseries.infdata.BW
    dmdelay = psr_utils.delay_from_DM(timeseries.infdata.DM, midfreq) - \
              psr_utils.delay_from_DM(timeseries.infdata.DM, hifreq)
    dmdelay_mjd = dmdelay/float(psr_utils.SECPERDAY)
    if debug:
        colour.cprint("High frequency (MHz): %f" % hifreq, 'debug')
        colour.cprint("Mid frequency (MHz): %f" % midfreq, 'debug')
        colour.cprint("DM delay added to TOAs (s): %g" % dmdelay, 'debug')
        colour.cprint("DM delay added to TOAs (MJD): %g" % dmdelay_mjd, 'debug')

    t0f = mjdf - phs*period/psr_utils.SECPERDAY + dmdelay_mjd
    t0i = mjdi
    shift,eshift,snr,esnr,b,errb,ngood,tphs = measure_phase(summed_pulse.profile, \
                            template_profile)
    # tphs is amount template is rotated by. It is originally 
    # measured in radians, convert to rotational phase
    tphs = tphs/(np.pi*2.0) % 1.0
    # tau and tau_err are the predicted phase of the pulse arrival
    tau, tau_err = shift/summed_pulse.N, eshift/summed_pulse.N
    if debug:
        colour.cprint("FFTFIT: Shift (bins): %f, Tau (phase): %f" % (shift, tau), 'debug')
        colour.cprint("FFTFIT: Shift error (bins): %f, Tau error (phase): %f" % \
                            (eshift, tau_err), 'debug')
    # Note: "error" flags are shift = 0.0 and eshift = 999.0
    if (np.fabs(shift) < 1e-7 and np.fabs(eshift-999.0) < 1e-7):
        raise FFTFitError("Error in FFTFIT. Bad return values.")
    # Send the TOA to STDOUT
    toaf = t0f + tau*period/float(psr_utils.SECPERDAY)
    if debug:
        colour.cprint("t0f (MJD): %r" % t0f, 'debug')
        colour.cprint("period (s): %r" % period, 'debug')
        colour.cprint("toaf (MJD): %r" % toaf, 'debug')
    newdays = int(np.floor(toaf))
    obs_code = telescopes.telescope_to_id[timeseries.infdata.telescope]
    psr_utils.write_princeton_toa(t0i+newdays, toaf-newdays, \
                                tau_err*period*1e6, midfreq, \
                                timeseries.infdata.DM, obs=obs_code)
    return tau, tphs


def measure_phase(profile, template):
    """Call FFTFIT on the profile and template to determine the
        following parameters: shift,eshift,snr,esnr,b,errb,ngood
        (returned as a tuple).  These are defined as in Taylor's
        talk at the Royal Society.

        pha1, the amount the template is rotated by (in radians)
        is also returned, in addition to the values mentioned
        above. pha1 is the last element in the tuple.
    """
    # This code is taken from Scott Ransom's PRESTO's get_TOAs.py
    c,amp,pha = fftfit.cprof(template)
    pha1 = pha[0]
    # Rotate the template
    pha = np.fmod(pha-np.arange(1,len(pha)+1)*pha1, 2.0*np.pi)
    shift,eshift,snr,esnr,b,errb,ngood = fftfit.fftfit(profile,amp,pha)
    return shift,eshift,snr,esnr,b,errb,ngood,pha1


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


def print_report(pulses, numpulses, nummasked, snrs=None, notes=None, \
                    quiet=False):
    """Print a report given the pulses provided.
    """
    print "Autopsy report:"
    print "\tTotal number of pulses searched: %s" % numpulses
    print "\tNumber of pulses thrown out: %s (%5.2f%%)" % (nummasked, \
                float(nummasked)/numpulses*100)
    print "\tNumber of good pulses found: %s (%5.2f%%)" % (len(pulses), \
                float(len(pulses))/numpulses*100)
    if len(pulses) > 0 and not quiet:
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
        pulse.plot(os.path.split(timeseries.basefn)[1], downfactor)


def joy_division_plot(pulses, timeseries, downfactor=1, hgt_mult=1):
    """Plot each pulse profile on the same plot separated
        slightly on the vertical axis.
        'timeseries' is the Datfile object dissected.
        Downsample profiles by factor 'downfactor' before plotting.
        hgt_mult is a factor to stretch the height of the paper.
    """
    first = True
    ppgplot.pgbeg("%s.joydiv.ps/CPS" % \
                    os.path.split(timeseries.basefn)[1], 1, 1)
    ppgplot.pgpap(10.25, hgt_mult*8.5/11.0) # Letter landscape
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
        # Interpolate before downsampling
        interp = ((copy_of_pulse.N/downfactor)+1)*downfactor
        copy_of_pulse.interpolate(interp)
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
        # Interpolate before downsampling
        interp = ((copy_of_pulse.N/downfactor)+1)*downfactor
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


def parse_on_pulse_regions(option, opt_str, value, parser):
    """Parse list of on-pulse regions from command line.
    """
    on_pulse = []
    for pair in value.split(','):
        lo, hi = pair.split(':')
        on_pulse.append((float(lo), float(hi)))
    setattr(parser.values, option.dest, on_pulse)

class FFTFitError(Exception):
    pass

if __name__ == '__main__':
    parser = optparse.OptionParser(usage="%prog --use-parfile PARFILE | --use-polycos POLYCOFILE | -p --period PERIOD [options] infile.dat", description="Given a input timeseries (a PRESTO .dat file) dissect it into individual pulses and record the pulses that surpass the signification threshold. Written by Patrick Lazarus.", version="%prog v0.9 (by Patrick Lazarus)", prog="dissect.py")
    parser.add_option('-t', '--threshold', dest='threshold', type='float', action='store', help="Only record pulses more significant than this threshold. (Default: 5).", default=5)
    parser.add_option('-n', '--no-output-files', dest='create_output_files', action='store_false', help="Do not create any output file for each significant pulse detected. (Default: create output files).", default=True)
    parser.add_option('--debug', dest='debug', action='store_true', help="Display debugging information. (Default: don't display debugging information).", default=False)
    parser.add_option('--no-text-files', dest='create_text_files', action='store_false', help="Do not create text file for each significant pulse detected. (Default: create text files).", default=True)
    parser.add_option('-q', '--quiet', dest='quiet', action='store_true', help="Output less information. (Default: Output full information).", default=False)
    parser.add_option('--no-toss', dest='no_toss', action='store_true', help="Do not toss out any pulses. (Default: Toss out partially masked profiles).", default=False)
    parser.add_option('-r', '--on-pulse-regions', dest='on_pulse_regions', type='string', action='callback', callback=parse_on_pulse_regions, help="Define (multiple) on-pulse regions. Beginning and end of each region should be separated by ':' and multiple pairs should be separated by ','. No spaces! Values should be given in terms of rotational phase, floats between 0.0 and 1.0. The on-pulse region is applied after the start of the observation has been shifted. (Default: None).", default=None)
    parser.add_option('-w', '--widths', dest='widths', type='string', action='callback', callback=parse_boxcar_widths, help="Boxcar widths (in number of samples) to use for smoothing profiles when searching for pulses. widths should be comma-separated _without_ spaces. (Default: Smooth with boxcar widths %s)" % DEFAULT_WIDTHS, default=DEFAULT_WIDTHS)
    parser.add_option('-s', '--shift-phase', dest='shift_phase', type='float', help="Set provided phase as the beginning of each pulse period. This is done by removing a piece of data from the beginning of the observation. If polycos, or parfile are used, then the phase is given by the polycos. If a constant period is used then the beginning of the observation is assumed to be phase=0.0. (Default: First pulse period begins at start of observation).", default=0.0)

    toa_group = optparse.OptionGroup(parser, "TOA Generation", "The following options are used to determine if/how TOAs are generated.")
    toa_group.add_option('--toas', dest='write_toas', action='store_true', help="Write TOAs to stdout. A TOA for each pulse will be written out unless --toa-threshold is provided, in which case consecutive pulses will be summed until sum profile's SNR surpases threshold.", default=False)
    toa_group.add_option('--template', dest='template', type='string', help="Required option if generating TOAs. This is the template profile to use.", default=None)
    toa_group.add_option('--toa-threshold', dest="toa_threshold", type='float', action='store', help="Threshold SNR for writing out TOAs. (Default: 0).", default=0)
    toa_group.add_option('--min-pulses', dest="min_pulses", type='int', action='store', help="Minimum number of pulses that must be added for writing out a single TOA. (Default: 1).", default=1)
    toa_group.add_option('--write-toa-files', dest='write_toa_files', action='store_true', help="Write out profiles used for TOAs as text files and postscript plots. (Default: Don't write out TOA profiles).", default=False)
    parser.add_option_group(toa_group)

    period_group = optparse.OptionGroup(parser, "Period Determination", "The following options are different methods for determine the spin period of the pulsar. Exactly one of these options must be provided.")
    period_group.add_option('--use-parfile', dest='parfile', type='string', action='store', help="Determine spin period from polycos generated by tempo using provided parfile.", default=None)
    period_group.add_option('--use-polycos', dest='polycofile', type='string', action='store', help="Determine spin period from polycos in the provided file.", default=None)
    period_group.add_option('-p', '--use-period', dest='period', type='float', action='store', help="Use constant period provided (in seconds).", default=None)
    parser.add_option_group(period_group)

    plot_group = optparse.OptionGroup(parser, "Plotting Options", "The following options affect only the output plots, not the data searching.")
    plot_group.add_option('-d', '--downsample', dest='downfactor', type='int', action='store', help="Down sample profiles by this factor before plotting. (Default: Don't downsample).", default=1)
    plot_group.add_option('--stretch-height', dest='heightstretch', type='float', action='store', help="Stretch height of JoyDiv plot by this factor. (Default: Do not stretch JoyDiv plot).", default=1)
    parser.add_option_group(plot_group)
    
    parser.add_option('--no-pulse-plots', dest='create_plot_files', action='store_false', help="Do not create plots for each significant pulse detected. (Default: create plots).", default=True)
    parser.add_option('--no-joydiv-plot', dest='create_joydiv_plot', action='store_false', help="Do not create Joy Division plot, where every profile is plotted in a single axes separated slightly in the vertical direction. (Default: create JoyDiv plot).", default=True)


    options, args = parser.parse_args()

    # Count number of period determination options are provided
    if (options.parfile is not None) + \
        (options.polycofile is not None) + \
        (options.period is not None) != 1:
        parser.print_usage()
        sys.stderr.write("Exactly one (1) period determination option must be provided! Exiting...\n\n")
        sys.exit(1)
    main()
