#!/usr/bin/env python

# A simple command line version of plotres written in python
# using matplotlib and numpy
#
#           Patrick Lazarus, Feb 26th, 2009

import pylab as p
import numpy as np
import residuals
import parfile as par
import binary_psr
import optparse
import sys
import slalib
import re
import os

class resid:
    def __init__(self):
        # Read TEMPO results (resid2.tmp, tempo.lis, timfile and parfiles)
        
        # Read residuals
        # Need to check if on 32-bit or 64-bit computer
        # (the following is a hack to find out if we're on a borg or borgii node)
        if os.uname()[4]=='i686':
            # 32-bit computer
            self.residuals = residuals.read_residuals()
        else:
            # 64-bit computer
            self.residuals = residuals.read_residuals_64bit()
        
        # Open tempo.lis. Parse it and find input .tim and .par files. Also find output .par file.
        inputfiles_re = re.compile(r"Input data from (.*\.tim),  Parameters from (.*\.par)")
        outputfile_re = re.compile(r"Assumed parameters -- PSR (.*)$")
        tempolisfile = open("tempo.lis")
        intimfn, inparfn, outparfn = None, None, None
        for line in tempolisfile:
            match = inputfiles_re.search(line)
            if match:
                intimfn = match.group(1).strip()
                inparfn = match.group(2).strip()
            else:
                match = outputfile_re.search(line)
                if match:
                    outparfn = "%s.par" % match.group(1).strip()
            if (intimfn != None) and (inparfn != None) and (outparfn != None):
                # Found what we're looking for no need to continue parsing the file
                break
        tempolisfile.close()

        # Read parfiles
        self.inpar = par.psr_par(inparfn)
        self.outpar = par.psr_par(outparfn)
        
        self.parse_options()
        
    def plot(self):
        xopt = self.options.xaxis.lower()
        if xopt == 'numtoa':
            xaxis = np.arange(self.residuals.numTOAs)
            xlabel = "TOA Number"
        elif xopt == 'mjd':
            xaxis = self.residuals.bary_TOA
            xlabel = "MJD"
        elif xopt == 'orbitphase':
            xaxis = self.residuals.orbit_phs
            xlabel = "Orbital Phase"
        elif xopt == 'year':
            xaxis = mjd_to_year(self.residuals.bary_TOA)
            xlabel = "Year"
        else:
            raise "How did you squeek through: options.xaxis = %s" % options.xaxis

        yopt = self.options.yaxis.lower()
        if self.options.postfit:
            title = "Postfit Residuals"
            if yopt == 'phase':
                yaxis = self.residuals.postfit_phs
                yerror = self.residuals.uncertainty/self.outpar.P0 # NOTE: Should use P at TOA not at PEPOCH
                ylabel = "Residuals (Phase)"
            elif yopt == 'usec':
                yaxis = self.residuals.postfit_sec*1e6
                yerror = self.residuals.uncertainty*1e6
                ylabel = "Residuals (uSeconds)"
            elif yopt == 'sec':
                yaxis = self.residuals.postfit_sec
                yerror = self.residuals.uncertainty
                ylabel = "Residuals (Seconds)"
            else:
                raise "How did you squeek through: options.yaxis = %s" % options.yaxis
        else:
            title = "Prefit Residuals"
            if yopt=='phase':
                yaxis = self.residuals.prefit_phs
                yerror = self.residuals.uncertainty/self.inpar.P0 # NOTE: Should use P at TOA not at PEPOCH
                ylabel = "Residuals (Phase)"
            elif yopt=='usec':
                yaxis = self.residuals.prefit_sec*1e6
                yerror = self.residuals.uncertainty*1e6
                ylabel = "Residuals (uSeconds)"
            elif yopt=='sec':
                yaxis = self.residuals.prefit_sec
                yerror = self.residuals.uncertainty
                ylabel = "Residuals (Seconds)"
            else:
                raise "How did you squeek through: options.yaxis = %s" % options.yaxis
        
        # Set up the plot
        fig = p.figure(figsize=(11,8.5))
        # Plot the residuals
        TOAcount = 0
        for lo,hi in self.freqbands:
            indices = (self.residuals.bary_freq>=lo) & (self.residuals.bary_freq<hi)
            p.errorbar(xaxis[indices], yaxis[indices], yerr=yerror[indices], fmt='.', label="%s - %s MHz" % (lo, hi))
            TOAcount += xaxis[indices].size
        # Finish off the plot
        p.axhline(0, ls='--', label="_nolegend_", c='k', lw=0.5)
        if self.options.mark_peri:
            if self.options.postfit:
                binpsr = binary_psr.binary_psr(self.outpar.FILE)
            else:
                binpsr = binary_psr.binary_psr(self.inpar.FILE)
            xmin, xmax = p.xlim()
            mjd_min = self.residuals.bary_TOA.min()
            mjd_max = self.residuals.bary_TOA.max()
            guess_mjds = np.arange(mjd_max + binpsr.par.PB, mjd_min - binpsr.par.PB, -binpsr.par.PB)
            for mjd in guess_mjds:
                peri_mjd = binpsr.most_recent_peri(float(mjd))
                if xopt == 'mjd':
                    p.axvline(peri_mjd, ls=':', label='_nolegend_', c='k', lw=0.5)
                elif xopt == 'year':
                    print "plotting peri passage"
                    p.axvline(mjd_to_year(peri_mjd), ls=':', label='_nolegend_', c='k', lw=0.5)
            p.xlim((xmin, xmax))
        p.xlabel(xlabel)
        p.ylabel(ylabel)
        p.title(title + " (Number of TOAs: %d)" % TOAcount)
        # Register event callback function and show the plot
        if self.options.legend:
            p.legend(loc='best')
        cid = fig.canvas.mpl_connect('key_press_event', keypress)
        if self.options.interactive:
            p.ion()
            p.show()
        else:
            # Save figure and quit
            savefigure()
            quit()

    def parse_options(self):
        if __name__ == '__main__':
            (self.options, arguments) = parser.parse_args()
        else:
            (self.options, arguments) = parser.parse_args([])
        
        if not self.options.freqs:
            freqbands = [['0', 'inf']]
        else:
            freqbands = []
            for fopt in self.options.freqs:
                f = fopt.split(':')
                if f[0]=='':
                    f[0] = '0'
                if f[-1]=='':
                    f[-1] = 'inf'
                if len(f) > 2:
                    for i in range(0, len(f)-1):
                        freqbands.append(f[i:i+2])
                else:
                    freqbands.append(f)
        freqbands = np.array(freqbands).astype(float)
        freqbands[freqbands.argsort(axis=0).transpose()[0]]
        if np.any(freqbands.flat != sorted(freqbands.flat)):
            raise "Frequency bands provided have overlaps or are inverted. Exiting!"
        self.freqbands = freqbands
       
        if not hasattr(self.outpar, 'BINARY'):
            # Cannot mark passage of periastron if pulsar is not in a binary
            self.options.mark_peri = False
       
        if self.options.xaxis.lower() not in ['numtoa', 'mjd', 'orbitphase', 'year']:
            raise "Option to -x/--x-axis (%s) is not permitted. Exiting!" % self.options.xaxis
        if self.options.yaxis.lower() not in ['phase', 'usec', 'sec']:
            raise "Option to -y/--y-axis (%s) is not permitted. Exiting!" % self.options.yaxis

def savefigure(savefn='./resid2.tmp.ps'):
        print "Saving plot to %s" % savefn
        p.savefig(savefn, orientation='landscape', papertype='letter')

def quit():
    print "Quiting..."
    sys.exit(0)

def keypress(event):
    if event.key.lower()=='q':
        quit()
    elif event.key.lower()=='s':
        savefigure()

def mjd_to_year(mjds):
    mjds = np.asarray(mjds)
    old_shape = mjds.shape # Remember original shape
    mjds.shape = (mjds.size, 1)
    years, months, days, fracs, stats = np.apply_along_axis(slalib.sla_djcl, 1, mjds).transpose()
    # Take into account leap years
    daysperyear = (((years % 4) == 0) & (((years % 100) != 0) | ((years % 400) == 0))) * 1 + 365.0
    years, days, stats = np.array([slalib.sla_clyd(*ymd) for ymd in np.vstack((years, months, days)).transpose()]).transpose()
    mjds.shape = old_shape # Change back to original shape
    return (years + (days + fracs) / daysperyear)

def main():
    r = resid()
    r.plot()
    
if __name__=='__main__':
    parser = optparse.OptionParser()
    parser.add_option('-f', '--freq', dest='freqs', action='append', help="Band of frequencies, in MHz, to be plotted (format xxx:yyy). Each band will have a different colour. Multiple -f/--freq options are allowed. (Default: Plot all frequencies in single colour.)", default=[])
    parser.add_option('-x', '--x-axis', dest='xaxis', type='string', help="Values to plot on x-axis. Must be one of {'numTOA', 'MJD', 'orbitphase', 'year'}. (Default: 'MJD')", default='MJD')
    parser.add_option('-y', '--y-axis', dest='yaxis', type='string', help="Values to plot on y-axis. Must be one of {'phase', 'usec', 'sec'}. (Default: 'phase')", default='phase')
    parser.add_option('--post', dest='postfit', action='store_true', help="Show postfit residuals. (Default: Plot postfit)", default=True)
    parser.add_option('--pre', dest='postfit', action='store_false', help="Show prefit residuals. (Default: Plot postfit.)")
    parser.add_option('-l', '--legend', dest='legend', action='store_true', help="Show legend of frequencies. (Default: Do not show legend.)", default=False)
    parser.add_option('--mark-peri', dest='mark_peri', action='store_true', help="Mark passage of periastron. (Default: don't mark periastron.)", default=False)
    parser.add_option('--non-interactive', dest='interactive', action='store_false', help="Save figure and exit. (Default: Show plot, only save if requested.)", default=True)
    
    main()


