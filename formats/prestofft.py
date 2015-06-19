"""
prestofft.py

Read PRESTO .fft file.
"""

import sys
import warnings
import os.path
import argparse
import types

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.interpolate

import infodata
import psr_utils

"""
Define object.
frequencies list (generator?)
powers list
normalisation method? Keep copy of powers as well as raw powers?
inf file
plot method
"""

COLOURS = ['r', 'b', 'g', 'm', 'c', 'y']

class PrestoFFT:
    def __init__(self, fftfn, inffn=None, maxfreq=None):
        """PrestoFFT object creator
            'fftfn' is filename of .fft file
            'inffn' is filename of .inf file
                If None the inffn will default to fftfn with .fft extension replaced by .inf
        """
        if not fftfn.endswith(".fft"):
            ValueError("FFT filename must end with '.fft'!")
        if not os.path.isfile(fftfn):
            ValueError("FFT file does not exist!\n\t(%s)" % fftfn)
        else:
            self.fftfn = fftfn
            self.fftfile = open(self.fftfn, "rb")

            if inffn is None:
                inffn = "%s.inf" % fftfn[:-4]
            else:
                inffn = inffn
            if not os.path.isfile(inffn):
                ValueError("Info file does not exist!\n\t(%s)" % inffn)

            self.inffn = inffn
            self.inf = infodata.infodata(inffn)

            freqs = np.fft.fftfreq(self.inf.N, self.inf.dt)
            self.freqs = freqs[freqs>=0]
            if maxfreq is not None:
                ntoread = np.sum(self.freqs<maxfreq)
                self.inf.N = ntoread
                self.freqs = self.freqs[:ntoread]
            else:
                ntoread = -1
            self.fft = self.read_fft(count=ntoread)
            self.phases = np.angle(self.fft)

            self.normalisation = "raw"
            self.powers = np.abs(self.fft)**2
            self.errs = None

    def close(self):
        self.fftfile.close()

    def interpolate(self, r, m=32):
        """Interpolate the value of the FFT at real bin indices 'r'.
            Use 'm' nearest bins when interpolating.

            Inputs:
                r: real bin indices to interpolate the FFT at.
                m: The number of nearby bins to use. Must be an
                    even integer.

            Output:
                interpfft: FFT coefficients interpolated at 'r'.
        """
        if (m % 2) is not 0:
            raise ValueError("Input 'm' must be an even integer: %s" % str(m))
        round_r = np.round(r).astype('int')
        k = round_r[:,np.newaxis]+np.arange(-m/2, m/2+1)
        coefs = self.fft[k]
        expterm = np.exp(-1.0j*np.pi*(r[:,np.newaxis]-k))
        sincterm = np.sinc(np.pi*(r[:,np.newaxis]-k))
        interpfft = np.sum(coefs*expterm*sincterm, axis=1)
        return interpfft

    def harmonic_sum(self, nharm=8):
        """Return the harmonically summed powers.

            Inputs:
                nharm: The number of harmonics to sum.

            Outputs:
                harmsummed: The harmonically summed powers.
        """
        nn = self.fft.size
        harmsummed = np.copy(self.powers[:nn/nharm])

        for nh in xrange(2, nharm+1):
            harmsummed += np.reshape(self.powers[:nn/nh*nh], (-1, nh))[:,0][:nn/nharm]
        return harmsummed

    def incoherent_harmonic_sum(self, nharm=8):
        """Return the harmonically summed FFT.

            Inputs:
                nharm: The number of harmonics to sum.

            Outputs:
                harmsummed: The harmonically summed FFT.
                freqs: The frequencies of each bin.
        """
        nn = self.fft.size
        harmsummed = np.copy(self.powers)
        
        for nh in xrange(2, 1+nharm):
            r = np.arange(nn)/float(nh)
            harmsummed += np.abs(self.interpolate(r, 2))**2
        return harmsummed, self.freqs/float(nharm)

    def coherent_harmonic_sum(self, nharm=8):
        """Return the harmonically summed FFT.

            Inputs:
                nharm: The number of harmonics to sum.

            Outputs:
                harmsummed: The harmonically summed FFT.
                freqs: The frequencies of each bin.
        """
        nn = self.fft.size
        harmsummed = np.copy(self.fft)
        
        for nh in xrange(2, 1+nharm):
            r = np.arange(nn)/float(nh)
            harmsummed += self.interpolate(r, 2)
        return np.abs(harmsummed)**2, self.freqs/float(nharm)

    def deredden(self, initialbuflen=6, maxbuflen=200):
        # Translated code from PRESTO's 'accel_utils.py'
        dered = np.copy(self.fft)
        # Takes care of the DC term
        dered[0] = 1 + 0j

        # Some initial values
        newoffset = 1
        fixedoffset = 1

        # Calculate initial values
        mean_old = np.median(self.powers[newoffset:newoffset+initialbuflen])/np.log(2.0)
        newoffset += initialbuflen
        lastbuflen = initialbuflen
        newbuflen = int(initialbuflen * np.log(newoffset))
        if newoffset > maxbuflen:
            newbuflen = maxbuflen
        
        while (newoffset + newbuflen) < len(dered):
            # Calculate the next mean
            mean_new = np.median(self.powers[newoffset:newoffset+newbuflen])/np.log(2.0)
            slope = (mean_new - mean_old) / (newbuflen + lastbuflen)

            # Correct the previous segment
            ioffs = np.arange(lastbuflen) # Index offsets for this block
            lineoffset = 0.5 * (newbuflen + lastbuflen)
            lineval = mean_old + slope * (lineoffset - ioffs)
            scaleval = 1.0/np.sqrt(lineval)
            dered[fixedoffset+ioffs] *= scaleval # multiplication is done 
                                                    # to both real and imag parts

            # Update our values
            fixedoffset += lastbuflen
            lastbuflen = newbuflen
            mean_old = mean_new
            newoffset += lastbuflen
            newbuflen = int(initialbuflen * np.log(newoffset))
            if newbuflen > maxbuflen:
                newbuflen = maxbuflen

        # Scale the last (partial) chunk the same way as the last point
        dered[fixedoffset:] *= scaleval[-1]
        return dered

    def estimate_power_errors(self, initialbuflen=6, maxbuflen=200, force=False):
        if not force and (self.errs is not None):
            return
        # Similar to deredden
        self.errs = np.zeros(len(self.fft))

        # Some initial values
        newoffset = 1
        fixedoffset = 1

        # Calculate initial values
        rms_old = np.std(self.powers[newoffset:newoffset+initialbuflen])
        newoffset += initialbuflen
        lastbuflen = initialbuflen
        newbuflen = int(initialbuflen * np.log(newoffset))
        if newoffset > maxbuflen:
            newbuflen = maxbuflen
        
        while (newoffset + newbuflen) < len(self.errs):
            # Calculate the next mean
            rms_new = np.std(self.powers[newoffset:newoffset+newbuflen])
            slope = (rms_new - rms_old) / (newbuflen + lastbuflen)

            # Correct the previous segment
            ioffs = np.arange(lastbuflen) # Index offsets for this block
            lineoffset = 0.5 * (newbuflen + lastbuflen)
            lineval = rms_old + slope * (lineoffset - ioffs)
            self.errs[fixedoffset+ioffs] =  lineval 

            # Update our values
            fixedoffset += lastbuflen
            lastbuflen = newbuflen
            rms_old = rms_new
            newoffset += lastbuflen
            newbuflen = int(initialbuflen * np.log(newoffset))
            if newbuflen > maxbuflen:
                newbuflen = maxbuflen

        # Assign the last (partial) chunk the same error as the last point
        self.errs[fixedoffset:] = lineval[-1]

    def fit_powers(self, freqlim=None, use_errors=True, **kwargs):
        # Use iminuit to fit power-law + DC to powers
        import iminuit
        if freqlim is None:
            freqlim = np.inf
            if self.inf.DM > 0:
                tdm = psr_utils.dm_smear(self.inf.DM, self.inf.BW, self.inf.lofreq+0.5*self.inf.BW)
                freqlim = 1.0/tdm
                print "Dispersion smearing time: %.2f ms" % (1000.0*tdm)
            freqlim = min(10.0, freqlim)
            print "Only fitting using frequencies up to %.2f Hz" % freqlim
        iuse = (self.freqs<freqlim)
        iuse[0] = 0 # Always ignore zeroth element
        
        if use_errors:
            # Compute power errors
            self.estimate_power_errors()

        # Define function to minimize
        def to_minimize(amp, index, dc):
            model = power_law(self.freqs[iuse], amp, index, dc) 
            diff = model-self.powers[iuse]
            #plt.figure()
            #plt.plot(self.freqs[1:], self.powers[1:], 'r-', alpha=0.5)
            #plt.plot(self.freqs[iuse], self.powers[iuse], 'k-', alpha=0.5)
            #plt.plot(self.freqs[iuse], model, 'k--', lw=2)
            #plt.xscale('log')
            #plt.yscale('log')
            #plt.xlabel("Freq (Hz)")
            #plt.ylabel("Raw Power")
            #plt.show()
            if use_errors:
                diff /= self.errs[iuse]
            return np.sum(diff**2)

        white = self.estimate_white_power_level(1000)
        kwargs.setdefault('amp', 1e14)
        kwargs.setdefault('error_amp', 1e12)
        kwargs.setdefault('limit_amp', (0, 1e20))
        kwargs.setdefault('index', -1.5)
        kwargs.setdefault('error_index', .1)
        kwargs.setdefault('limit_index', (-10.0, 0.0))
        kwargs.setdefault('dc', white)
        kwargs.setdefault('error_dc', white*0.1)
        #kwargs.setdefault('limit_dc', (white*0.01, np.max(self.powers[1:])))
        kwargs.setdefault('fix_dc', True)
        m = iminuit.Minuit(to_minimize, 
                           frontend=iminuit.ConsoleFrontend,
                           print_level=0, **kwargs)
        # Minimize
        m.migrad()
        return m.values

    def estimate_white_power_level(self, minfreq=1000):
        """Estimate the white noise level in the power spectrum.
            This is the median of the powers above a given frequency.

            Inputs:
                minfreq: Lower frequency limit for 'white' interval.
                    (Default: 1000 Hz)

            Outputs:
                white: The white noise level of the power spectrum.
        """
        return np.median(self.powers[self.freqs>minfreq])

    def plot_power_fit(self, powerlaws):
        #plt.plot(self.freqs[1:], self.powers[1:], 'k-')
        for ii, (amp, index, dc) in enumerate(powerlaws):
            c = COLOURS[ii%len(COLOURS)]
            model = power_law(self.freqs, amp, index, dc)
            plt.plot(self.freqs[1:], model[1:], ls='--', c=c,
                     label=r"A=%.2g, $\alpha$=%.3g, DC=%.2g" % (amp, index,dc))
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power")
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(loc='upper right', prop=dict(size='x-small'))

    def read_fft(self, count=-1):
        """Read 'count' powers from .fft file and return them.
            power = real*real + imag*imag
        """
        fft = np.fromfile(self.fftfile, dtype=np.dtype('c8'), count=count)
        return fft

    def plot(self, **kwargs):
        """Plot the power spectrum.
        """
        plt.plot(self.freqs, self.powers, **kwargs)
        plt.title(self.fftfn)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power")
        plt.xscale('log')
        plt.yscale('log')

    def plot_3pane(self):
        """Plot the power spectrum in 3 panes.
            If 'showzap' is True show PALFA zaplist.
        """
        ones = (self.freqs>=1) & (self.freqs<10)
        tens = (self.freqs>=10) & (self.freqs<100) 
        hundreds = (self.freqs>=100) & (self.freqs<1000)

        fig = plt.figure(figsize=(10,8))
        plt.subplots_adjust(hspace=0.25)
    
        # plot mean power spectrum (in three parts)
        axones = plt.subplot(3,1,1) 
        plt.plot(self.freqs[ones], self.powers[ones], 'k-', lw=0.5)
        plt.ylabel("Power") 
        plt.xscale('log')

        axtens = plt.subplot(3,1,2, sharey=axones)
        plt.plot(self.freqs[tens], self.powers[tens], 'k-', lw=0.5)
        plt.ylabel("Power")
        plt.xscale('log')

        axhundreds = plt.subplot(3,1,3, sharey=axones)
        plt.plot(self.freqs[hundreds], self.powers[hundreds], 'k-', lw=0.5)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power")
        plt.xscale('log')

        maxpwr = np.max(self.powers[(self.freqs>=1) & (self.freqs<1000)])
        axones.set_ylim(0, maxpwr*1.1)

        plt.suptitle("Power Spectrum (%s)" % self.fftfn)

    def plot_zaplist(self, zapfile, fc='b', ec='none', 
                     alpha=0.25, zorder=-1, **kwargs):
        # Plot regions that are zapped
        zaplist = np.loadtxt(zapfile)
        ax = plt.gca()
        for freq, width in zaplist:
            plt.axvspan(freq-width/2.0, freq+width/2.0, \
                        fill=True, fc=fc, ec=ec, \
                        alpha=alpha, zorder=zorder, **kwargs)
        plt.figtext(0.025, 0.03, "Zaplist file: %s" % zapfile, size="xx-small")


def power_law(freqs, amp, index, dc):
    rednoise = amp*freqs**index
    return rednoise + dc


def get_smear_response(ddm, **obs):
    if ddm != 0:
        bw = obs['chan_width']*obs['numchan']
        fhi = obs['lofreq']+bw
        smear = smearing_function(obs['lofreq'], fhi, ddm, obs.get('bandpass', None))
        times = np.arange(obs['N'])*obs['dt']
        weights = smear(times)
        weights /= np.sum(weights)

        freqs = np.fft.fftfreq(obs['N'], obs['dt'])
        freqs = freqs[freqs>=0]
        fft = np.fft.rfft(weights)[:len(freqs)]
        response = scipy.interpolate.interp1d(freqs, (np.abs(fft))**2)
    else:
        response = lambda freq: 1
    return response
    

def smearing_function(flo, fhi, ddm, bandpass=None):
    """Return function that returns the kernel used to smear data due to a wrong DM

        Inputs:
            flo: Low frequency of observing band (MHz)
            fhi: High frequency of observing band (MHz)
            ddm: Difference in DM with respect to the DM of the signal (pc/cc)

        Output:
            smear: A function that returns the smearing kernel.
    """
    if bandpass is not None:
        freqs = np.linspace(flo,fhi, len(bandpass))
        delay = 4.15e3*ddm*(freqs**-2 - fhi**-2)
        isort = np.argsort(delay)
        bandpass[~np.isfinite(bandpass)] = 0
        interp = scipy.interpolate.interp1d(delay[isort], bandpass[isort], 
                                            bounds_error=False, fill_value=0)
    else:
        interp = lambda time: 1

    tmax = 4.15e3*ddm*(flo**-2 - fhi**-2)
    def smear(times):
        weights = interp(times)/np.sqrt(times/4.15e3/ddm + fhi**-2)/(2*4.15e3*ddm)
        if tmax > 0:
            weights[(times<0) | (tmax<times)] = 0
        else:
            weights[(times<tmax) | (0<times)] = 0
        return weights
    return smear


def main():
    pfftfn = args.fftfn
    pfft = PrestoFFT(pfftfn, maxfreq=args.max_freq)

    plparams = args.powerlaw
    if args.do_fit:
        powerlaw = pfft.fit_powers(freqlim=args.freq_lim)
        print "Amplitude:", powerlaw['amp']
        print "Index:", powerlaw['index']
        print "DC offset:", powerlaw['dc']
        plparams.append((powerlaw['amp'], powerlaw['index'], powerlaw['dc']))
        write_powerlaw_to_file(powerlaw, pfftfn, outfn=(pfftfn+".plp"))
    if args.do_plot:
        fig = plt.figure()
        pfft.plot(c='k', lw=0.5)
        plt.xscale(args.xscale)
        plt.yscale(args.yscale)
        if args.zapfile is not None:
            pfft.plot_zaplist(args.zapfile)
        if plparams:
            pfft.plot_power_fit(plparams)
        if args.save_plot:
            plt.savefig(pfftfn+".png") 
        if args.show_plot:
            def close(event):
                if event.key in ('q','Q'):
                    plt.close()
            fig.canvas.mpl_connect("key_press_event", close)
            plt.show()


def write_powerlaw_to_file(powerlaw, pfftfn, outfn):
    with open(outfn, 'w') as ff:
        ff.write("# Power law parameters for %s\n" % pfftfn)
        ff.write("# Amplitude   Index   DC offset\n")
        ff.write("%(amp)g\t%(index)g\t%(dc)g\n" % powerlaw)


class PowerLawFromFile(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        plparams = getattr(namespace, self.dest)
        with open(values, 'r') as ff:
            for line in ff:
                line = line.partition("#")[0].strip()
                if not line:
                    continue
                split = line.split()
                if len(split) != 3:
                    raise ValueError("Each line in power law file must "
                                     "contain 3 params (amp, index, dc). "
                                     "%d params in:\n%s" % (len(split), line))
                plparams.append(tuple([float(x) for x in split]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Work with PRESTO FFT files.")
    parser.add_argument("--fit", dest='do_fit', action='store_true', \
                        help="Include a power-law fit to the power spectrum "
                             "when plotting.")
    parser.add_argument("--no-plot", dest='do_plot', action='store_false', \
                        help="Make plot of power spectrum.")
    parser.add_argument("--save-plot", dest="save_plot", action='store_true', \
                        help="Save plot to file. The figure will "
                             "be saved to <fft fn>.png")
    parser.add_argument("--no-show-plot", dest="show_plot", action='store_false', \
                        help="Show plot of power spectrum.")
    parser.add_argument("--fit-freq-lim", dest='freq_lim', type=float, \
                        default=10.0,
                        help="When fitting the power spectrum's red noise, only "
                             "use frequencies less than this value. " 
                             "(Default: minimum of 10 Hz and 1/DM smearing time)")
    parser.add_argument("--powerlaw", dest='powerlaw', action='append', \
                        type=float, nargs=3, default=[], \
                        help="Parameters for power law to plot. Three "
                             "parameters must be provided: amplitude index DC")
    parser.add_argument("--pl-params", dest='powerlaw', action=PowerLawFromFile, \
                        type=str, default=[], \
                        help="File with parameters for power law to plot. Three "
                             "parameters must be provided: amplitude index DC")
    parser.add_argument("--xscale", type=str, dest='xscale',
                        default='log',
                        help="Matplotlib-style scaling type for the X-axis. "
                             "(Default: log)")
    parser.add_argument("--yscale", type=str, dest='yscale',
                        default='log',
                        help="Matplotlib-style scaling type for the Y-axis. "
                             "(Default: log)")
    parser.add_argument("-z", "--zapfile", type=str, dest='zapfile', \
                        help="Zap file to show when plotting.")
    parser.add_argument("--max-freq", dest="max_freq", type=float,
                        help="Maximum frequency to plot. (Default: plot all)")
    parser.add_argument("fftfn", type=str,
                        help="PRESTO FFT file to use. The corresponding *.inf " \
                             "file must be present.")
    args = parser.parse_args()
    main()
