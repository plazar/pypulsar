import optparse
import sys
import warnings

import numpy as np
import matplotlib.pyplot as plt
import psr_utils
import fftfit
import pulse
import telescopes

def main():
    pulses = get_pulses(args)
    # Sort pulses by increasing MJD
    cmp_mjd = lambda this, other: cmp(this.mjd, other.mjd)
    pulses.sort(cmp=cmp_mjd)

    print "Number of pulses loaded: %d" % len(pulses)
    if len(pulses):
        numtoas = 0
        print "Generating TOAs. Please wait..."
        print "TOA threshold:", options.toa_threshold
        print "Min number of pulses for a TOA:", options.min_pulses
        print "Profile template used:", options.template
        # Extract profile template file
        # (Note: first column is index, second column is profile values)
        template = np.loadtxt(options.template, usecols=(1,))
        # Get TOAs and write them to stdout
        summed_pulse = None
        numpulses = 0
        registry = [] # too keep track of what pulses make up what TOA
        for pulse in pulses:
            if summed_pulse is None:
                summed_pulse = pulse.make_copy()
                numpulses = 1
                pulsenums = [pulse.number]
            else:
                summed_pulse= add_profs(summed_pulse, pulse)
                numpulses += 1
                pulsenums.append(pulse.number)
            if numpulses < options.min_pulses:
                continue
            if get_snr(summed_pulse) > options.toa_threshold:
                # Interpolate and downsample summed_pulse so
                # it is same size as template profile
                summed_pulse.interp_and_downsamp(template.size)
                summed_pulse.scale()
                pulseshift, templateshift = write_toa(summed_pulse, \
                            template, options.debug)
                numtoas += 1
                if options.write_toa_files:
                    # TOA profiles are already downfactored
                    # do not downfactor more when creating plot
                    plot_toa(numtoas, summed_pulse, template, \
                            pulseshift, templateshift)
                    summed_pulse.write_to_file("TOA%d" % numtoas)
                registry.append(pulsenums)
                summed_pulse = None
                numpulses = 0
        print "Number of TOAs: %d" % numtoas
        print "Number of pulses thrown out because 'min pulses' requirement " \
                "or SNR threshold not met: %d" % numpulses
        print "Registry of pulses:"
        for i, pulsenums in enumerate(registry):
            print "\tTOA #%d: %d pulses" % (i+1, len(pulsenums))
            print "\t\t(pulse numbers: %s)" % ', '.join([str(n) for n in pulsenums])


def get_pulses(pulsefiles):
    """Given a list of pulse files
        return a list of Pulse objects.
    """
    pulses = []
    for pf in pulsefiles:
        pulses.append(pulse.read_pulse_from_file(pf))
    return pulses


def add_profs(pulse, other):
    """Given two Pulse objects add the profile of 'other'
        to the profile of 'pulse'.

        Return a copy of 'pulse' with the summed profile.

        NOTE: The summed profile (of 'pulse') will have the
                length of the shorter of the two input profiles.
              No information of 'pulse' is changed.
    """
    copy = pulse.make_copy()
    # Truncate to size of smaller profile
    N = np.min([len(pulse.profile), len(other.profile)])
    copy.profile = pulse.profile[:N] + other.profile[:N]
    copy.N = N
    copy.duration = N*copy.dt
    return copy


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


def write_toa(summed_pulse, template_profile, debug=False):
    """Given a SummedPulse generate a TOA and write it to stdout. 
        'template_profile' is simply a numpy array.
        'debug' is a boolean value that determines if debugging
        info should be displayed.
        Returns shift required to line up template and pulse.
    """
    if template_profile is None:
        raise ValueError("A template profile MUST be provided.")
    # This code is taken from Scott Ransom's PRESTO's get_TOAs.py
    mjdi = int(summed_pulse.mjd) # integer part of MJD
    mjdf = summed_pulse.mjd - mjdi # fractional part of MJD
    # Period is duration of profile in seconds
    period = summed_pulse.dt*len(summed_pulse.profile)
    
    # Caclulate offset due to shifting channels to account for DM
    # Hifreq doesn't have a half-channel offset 
    # (why? because get_TOAs.py doesn't. Why...)
    # Why subtract 1 channel to get hifreq?
    hifreq = summed_pulse.lofreq + summed_pulse.bw - summed_pulse.chan_width
    midfreq = summed_pulse.lofreq - 0.5*summed_pulse.chan_width + 0.5*summed_pulse.bw
    dmdelay = psr_utils.delay_from_DM(summed_pulse.dm, midfreq) - \
              psr_utils.delay_from_DM(summed_pulse.dm, hifreq)
    dmdelay_mjd = dmdelay/float(psr_utils.SECPERDAY)
    if debug:
        colour.cprint("High frequency (MHz): %f" % hifreq, 'debug')
        colour.cprint("Mid frequency (MHz): %f" % midfreq, 'debug')
        colour.cprint("DM delay added to TOAs (s): %g" % dmdelay, 'debug')
        colour.cprint("DM delay added to TOAs (MJD): %g" % dmdelay_mjd, 'debug')

    t0f = mjdf + dmdelay_mjd
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
    obs_code = telescopes.telescope_to_id[summed_pulse.telescope]
    psr_utils.write_princeton_toa(t0i+newdays, toaf-newdays, \
                                tau_err*period*1e6, midfreq, \
                                summed_pulse.dm, obs=obs_code)
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


class FFTFitError(Exception):
    pass


if __name__ == '__main__':
    parser = optparse.OptionParser(usage="%prog [options] {PROFILE FILES}", \
                    description="Write TOAs to stdout. A TOA for each pulse " \
                         "will be written out unless --toa-threshold is " \
                         "provided, in which case consecutive pulses will " \
                         "be summed until the summed profile's SNR surpases " \
                         "the given threshold.", \
                    version="%prog v0.9", prog="pulses_to_toas.py")
    parser.add_option('--debug', dest='debug', action='store_true', \
                    help="Display debugging information. (Default: don't " \
                         "display debugging information).", \
                    default=False)
    toa_group = optparse.OptionGroup(parser, "TOA Generation", \
                                "The following options are used to determine if/how " \
                                "TOAs are generated.")
    toa_group.add_option('--template', dest='template', type='string', \
                    help="Required option if generating TOAs. This is the " \
                         "template profile to use.", default=None)
    toa_group.add_option('--toa-threshold', dest="toa_threshold", type='float',  \
                    help="Threshold SNR for writing out TOAs. (Default: 0).", \
                    action='store', default=0)
    toa_group.add_option('--min-pulses', dest="min_pulses", type='int', \
                    help="Minimum number of pulses that must be added for " \
                         "writing out a single TOA. (Default: 1).", \
                    action='store', default=1)
    toa_group.add_option('--write-toa-files', dest='write_toa_files', \
                    help="Write out profiles used for TOAs as text files and " \
                         "postscript plots. (Default: Don't write out TOA profiles).", \
                    action='store_true', default=False)
    parser.add_option_group(toa_group)

    options, args = parser.parse_args()
    main()
