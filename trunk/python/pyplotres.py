#!/usr/bin/env python

# A simple command line version of plotres written in python
# using matplotlib and numpy
#
#           Patrick Lazarus, Feb 26th, 2009

import optparse
import sys
import re
import os
import types
import warnings

import matplotlib.pyplot as plt
import numpy as np

import slalib
import binary_psr
import parfile as par
import residuals

class TempoResults:
    def __init__(self, freqbands=[[0, 'inf']]):
        """Read TEMPO results (resid2.tmp, tempo.lis, timfile and parfiles)
            freqbands is a list of frequency pairs to display.
        """
        # Open tempo.lis. Parse it and find input .tim and .par files. Also find output .par file.
        inputfiles_re = re.compile(r"Input data from (.*\.tim.*),  Parameters from (.*\.par.*)")
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

        # Record filename
        self.inparfn = inparfn
        self.outparfn = outparfn
        self.intimfn = intimfn
        
        # Read parfiles
        self.inpar = par.psr_par(inparfn)
        self.outpar = par.psr_par(outparfn)
        
        # Read residuals
        # Need to check if on 32-bit or 64-bit computer
        # (the following is a hack to find out if we're on a borg or borgii node)
        if os.uname()[4]=='i686':
            # 32-bit computer
            r = residuals.read_residuals()
        else:
            # 64-bit computer
            r = residuals.read_residuals_64bit()

        self.max_TOA = r.bary_TOA.max()
        self.min_TOA = r.bary_TOA.min()
        self.freqbands = freqbands 
        self.residuals = {}
        for lo,hi in self.freqbands:
            indices = (r.bary_freq>=lo) & (r.bary_freq<hi)
            self.residuals[get_freq_label(lo, hi)] = \
                 Resids(r.bary_TOA[indices], r.bary_freq[indices], \
                        np.arange(r.numTOAs)[indices], r.orbit_phs[indices], \
                        r.postfit_phs[indices], r.postfit_sec[indices], \
                        r.prefit_phs[indices], r.prefit_sec[indices], \
                        r.uncertainty[indices], r.weight[indices], \
                        self.inpar, self.outpar)


class Resids:
    """The Resids object contains the following information
        about TEMPO residuals:
            bary_TOA
            bary_freq
            numTOAs
            orbit_phs
            postfit_phs
            postfit_sec
            prefit_phs
            prefit_sec
            uncertainty
            weight
    """
    def __init__(self, bary_TOA, bary_freq, TOA_index, orbit_phs, \
                    postfit_phs, postfit_sec, prefit_phs, prefit_sec, \
                    uncertainty, weight, inpar, outpar):
        self.bary_TOA = bary_TOA
        self.bary_freq = bary_freq
        self.TOA_index = TOA_index
        self.orbit_phs = orbit_phs
        self.postfit_phs = postfit_phs
        self.postfit_sec = postfit_sec
        self.prefit_phs = prefit_phs
        self.prefit_sec = prefit_sec
        self.uncertainty = uncertainty
        self.weight = weight
        self.inpar = inpar
        self.outpar = outpar
        

    def get_xdata(self, key):
        """Return label describing xaxis and the corresponding 
            data given keyword 'key'.
        """
        if not isinstance(key, types.StringType):
            raise ValueError("key must be of type string.")
        xopt = key.lower()
        if xopt == 'numtoa':
            xdata = self.TOA_index
            xlabel = "TOA Number"
        elif xopt == 'mjd':
            xdata = self.bary_TOA
            xlabel = "MJD"
        elif xopt == 'orbitphase':
            xdata = self.orbit_phs
            xlabel = "Orbital Phase"
        elif xopt == 'year':
            xdata = mjd_to_year(self.bary_TOA)
            xlabel = "Year"
        else:
            raise ValueError("Unknown xaxis type (%s)." % xopt)
        return (xlabel, xdata)

    
    def get_ydata(self, key, postfit=True):
        """Return label describing yaxis and the corresponding 
            data/errors given keyword 'key'.
            'postfit' is a boolean argument that determines if
            postfit, or prefit data is to be returned.
        """
        if not isinstance(key, types.StringType):
            raise ValueError("key must be of type string.")
        yopt = key.lower()
        if postfit:
            if yopt == 'phase':
                ydata = self.postfit_phs
                #
                # NOTE: Should use P at TOA not at PEPOCH
                #
                yerror = self.uncertainty/self.outpar.P0 
                ylabel = "Residuals (Phase)"
            elif yopt == 'usec':
                ydata = self.postfit_sec*1e6
                yerror = self.uncertainty*1e6
                ylabel = "Residuals (uSeconds)"
            elif yopt == 'sec':
                ydata = self.postfit_sec
                yerror = self.uncertainty
                ylabel = "Residuals (Seconds)"
            else:
                raise ValueError("Unknown yaxis type (%s)." % yopt)
        else:
            if yopt=='phase':
                ydata = self.prefit_phs
                #
                # NOTE: Should use P at TOA not at PEPOCH
                #
                yerror = self.uncertainty/self.inpar.P0 
                ylabel = "Residuals (Phase)"
            elif yopt=='usec':
                ydata = self.prefit_sec*1e6
                yerror = self.uncertainty*1e6
                ylabel = "Residuals (uSeconds)"
            elif yopt=='sec':
                ydata = self.prefit_sec
                yerror = self.uncertainty
                ylabel = "Residuals (Seconds)"
            else:
                raise ValueError("Unknown yaxis type (%s)." % yopt)
        return (ylabel, ydata, yerror)


def plot_data(tempo_results, xkey, ykey, postfit=True, prefit=False, \
            interactive=True, mark_peri=False, show_legend=True):
    # figure out what should be plotted
    # True means to plot postfit
    # False means to plot prefit
    if postfit and prefit:
        to_plot_postfit = [False, True]
    elif postfit and not prefit:
        to_plot_postfit = [True]
    elif not postfit and prefit:
        to_plot_postfit = [False]
    else:
        raise ValueError("At least one of prefit and postfit must be True.")
    subplot = 1
    numsubplots = len(to_plot_postfit)
    global axes
    axes = []
    handles = []
    labels = []
    for usepostfit in to_plot_postfit:
        TOAcount = 0
        # All subplots are in a single column
        if subplot == 1:
            axes.append(plt.subplot(numsubplots, 1, subplot))
        else:
            axes.append(plt.subplot(numsubplots, 1, subplot, sharex=axes[0]))
        
        for lo,hi in tempo_results.freqbands:
            freq_label = get_freq_label(lo, hi)
            resids = tempo_results.residuals[freq_label]
            xlabel, xdata = resids.get_xdata(xkey)
            ylabel, ydata, yerr = resids.get_ydata(ykey, usepostfit)
            if len(xdata):
                # Plot the residuals
                handle = plt.errorbar(xdata, ydata, yerr=yerr, fmt='.', \
                                        label=freq_label, picker=5)
                if subplot == 1:
                    handles.append(handle[0])
                    labels.append(freq_label)
                TOAcount += xdata.size
        # Finish off the plot
        plt.axhline(0, ls='--', label="_nolegend_", c='k', lw=0.5)
        axes[-1].ticklabel_format(style='plain', axis='x')
       
        if mark_peri and hasattr(tempo_results.outpar, 'BINARY'):
            # Be sure to check if pulsar is in a binary
            # Cannot mark passage of periastron if not a binary 
            if usepostfit:
                binpsr = binary_psr.binary_psr(tempo_results.outpar.FILE)
            else:
                binpsr = binary_psr.binary_psr(tempo_results.inpar.FILE)
            xmin, xmax = plt.xlim()
            mjd_min = tempo_results.min_TOA
            mjd_max = tempo_results.max_TOA
            guess_mjds = np.arange(mjd_max + binpsr.par.PB, \
                                mjd_min - binpsr.par.PB, -binpsr.par.PB)
            for mjd in guess_mjds:
                peri_mjd = binpsr.most_recent_peri(float(mjd))
                if xopt == 'mjd':
                    plt.axvline(peri_mjd, ls=':', label='_nolegend_', c='k', lw=0.5)
                elif xopt == 'year':
                    print "plotting peri passage"
                    plt.axvline(mjd_to_year(peri_mjd), ls=':', label='_nolegend_', c='k', lw=0.5)
            plt.xlim((xmin, xmax))
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if usepostfit:
            plt.title("Postfit Redisuals (Number of TOAs: %d)" % TOAcount)
        else:
            plt.title("Prefit Redisuals (Number of TOAs: %d)" % TOAcount)
        subplot += 1
    
    if numsubplots > 1:
        # Increase spacing between subplots.
        plt.subplots_adjust(hspace=0.25)
   
    # Write name of input files used for timing on figure
    if interactive:
        fntext = "TOA file: %s, Parameter file: %s" % \
                    (tempo_results.intimfn, tempo_results.inparfn)
        figure_text = plt.figtext(0.01, 0.01, fntext, verticalalignment='bottom', \
                            horizontalalignment='left')
    
    if show_legend:
        leg = plt.figlegend(handles, labels, 'upper right')
        leg.legendPatch.set_alpha(0.5)
    

def create_plot(tempo_results, xkey, ykey, postfit=True, prefit=False, \
            interactive=True, mark_peri=False, show_legend=True):
    # Set up the plot
    fig = plt.figure(figsize=(11,8.5))
   
    plot_data(tempo_results, xkey, ykey, postfit, prefit, \
            interactive, mark_peri, show_legend)
    
    # Register event callbacks function and show the plot
    cid_keypress = fig.canvas.mpl_connect('key_press_event', keypress)
    cid_pick = fig.canvas.mpl_connect('pick_event', pick)
    if interactive:
        plt.ion()
        plt.show()
    else:
        # Save figure and quit
        savefigure()
        quit()


def get_freq_label(lo, hi):
    """Return frequency label given a lo and hi
        frequency pair.
    """
    return "%s - %s MHz" % (lo, hi)


def savefigure(savefn='./resid2.tmp.ps'):
        print "Saving plot to %s" % savefn
        plt.savefig(savefn, orientation='landscape', papertype='letter')


def quit():
    print "Quiting..."
    sys.exit(0)


def pick(event):
    global tempo_results
    index = event.ind
    if event.mouseevent.inaxes:
        ylabel = event.mouseevent.inaxes.get_ylabel()
        title = event.mouseevent.inaxes.get_title()
    if len(index) == 1:
        freq_label = event.artist.get_label()
        r = tempo_results.residuals[freq_label]
        print "TOA Selected:"
        print "\tNumber:", r.TOA_index[index][0]
        print "\tEpoch (MJD):", r.bary_TOA[index][0]
        if "(Phase)" in ylabel:
            print "\tPre-fit residual (phase):", r.prefit_phs[index][0]
            print "\tPost-fit residual (phase):", r.postfit_phs[index][0]
            if "Prefit" in title:
                print "\tUncertainty (phase):", r.uncertainty[index][0]/r.inpar.P0
            elif "Postfit" in title:
                print "\tUncertainty (phase):", r.uncertainty[index][0]/r.outpar.P0
            else:
                raise ValueError("Cannot determine pre/post-fit from title (%s)" % title)
        elif "(uSeconds)" in ylabel:
            print "\tPre-fit residual (usec):", r.prefit_sec[index][0]*1e6
            print "\tPost-fit residual (usec):", r.postfit_sec[index][0]*1e6
            print "\tUncertainty (usec):", r.uncertainty[index][0]*1e6
        elif "(Seconds)" in ylabel:
            print "\tPre-fit residual (sec):", r.prefit_sec[index][0]
            print "\tPost-fit residual (sec):", r.postfit_sec[index][0]
            print "\tUncertainty (sec):", r.uncertainty[index][0]
        else:
            raise ValueError("Unknown unit for y-axes (%s)" % ylabel)
        print "\tFrequency (MHz):", r.bary_freq[index][0]
    else:
        print "Multiple TOAs selected. Zoom in and try again."


def keypress(event):
    global tempo_results
    global options
    if type(event.key) == types.StringType:
        
        if event.key.lower() == 'q':
            quit()
        elif event.key.lower() == 's':
            savefigure()
        elif event.key.lower() == 'z':
            # Turn on zoom mode
            print "Toggling zoom mode..."
            event.canvas.toolbar.zoom()
        elif event.key.lower() == 'h':
            # Restore plot to original view
            print "Restoring plot..."
            event.canvas.toolbar.home()
        elif event.key.lower() == ',' or event.key.lower() == '<':
            # Go back to previous plot view
            print "Going back..."
            event.canvas.toolbar.back()
        elif event.key.lower() == '.' or event.key.lower() == '>':
            # Go forward to next plot view
            print "Going forward..."
            event.canvas.toolbar.forward()
        elif event.key.lower() == 'x':
            # Set x-axis limits
            print "Setting x-axis limits. User input required..."
            xmin = raw_input("X-axis minimum: ")
            xmax = raw_input("X-axis maximum: ")
            try:
                xmin = float(xmin)
                xmax = float(xmax)
                if xmax <= xmin:
                    raise ValueError
            except ValueError:
                print "Bad values provided!"
                return
            plt.xlim(xmin, xmax)
        elif event.key.lower() == 'y':
            global axes
            # Set y-axis limits
            print "Setting y-axis limits. User input required..."
            if len(axes) == 2:
                axes_to_adjust = raw_input("Axes to adjust (pre/post): ")
                if axes_to_adjust.lower().startswith('pre'):
                    plt.axes(axes[0])
                elif axes_to_adjust.lower().startswith('post'):
                    plt.axes(axes[1])
                else:
                    raise ValueError
            ymin = raw_input("Y-axis minimum: ")
            ymax = raw_input("Y-axis maximum: ")
            try:
                ymin = float(ymin)
                ymax = float(ymax)
                if ymax <= ymin:
                    raise ValueError
            except ValueError:
                print "Bad values provided!"
                return
            plt.ylim(ymin, ymax)
        elif event.key == ' ':
            # Reload residuals and replot
            print "Reloading..."
            plt.clf() # clear figure
            tempo_results = TempoResults(options.freqbands)
            plot_data(tempo_results, options.xaxis, options.yaxis, 
                    options.postfit, options.prefit, options.mark_peri, \
                    options.legend)
            

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


def parse_options():
    (options, sys.argv) = parser.parse_args()
    
    if not options.freqs:
        freqbands = [['0', 'inf']]
    else:
        freqbands = []
        for fopt in options.freqs:
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
        raise ValueError("Frequency bands have overlaps or are inverted.")
    options.freqbands = freqbands
   
    options.mark_peri = False
   
    if not options.prefit and not options.postfit:
        # If neither prefit or postfit are selected
        # show postfit
        options.postfit = True
   
    if options.xaxis.lower() not in ['numtoa', 'mjd', 'orbitphase', 'year']:
        raise BadOptionValueError("Option to -x/--x-axis (%s) is not permitted." % \
                            options.xaxis)
    if options.yaxis.lower() not in ['phase', 'usec', 'sec']:
        raise BadOptionValueError("Option to -y/--y-axis (%s) is not permitted." % \
                            options.yaxis)
    return options


def main():
    global tempo_results
    global options
    options = parse_options()
    tempo_results = TempoResults(options.freqbands)
    create_plot(tempo_results, options.xaxis, options.yaxis, options.postfit, \
            options.prefit, options.interactive, options.mark_peri, \
            options.legend)


class BadOptionValueError(ValueError):
    """Bad value passed to option parser.
    """
    pass


if __name__=='__main__':
    parser = optparse.OptionParser()
    parser.add_option('-f', '--freq', dest='freqs', action='append', help="Band of frequencies, in MHz, to be plotted (format xxx:yyy). Each band will have a different colour. Multiple -f/--freq options are allowed. (Default: Plot all frequencies in single colour.)", default=[])
    parser.add_option('-x', '--x-axis', dest='xaxis', type='string', help="Values to plot on x-axis. Must be one of {'numTOA', 'MJD', 'orbitphase', 'year'}. (Default: 'MJD')", default='MJD')
    parser.add_option('-y', '--y-axis', dest='yaxis', type='string', help="Values to plot on y-axis. Must be one of {'phase', 'usec', 'sec'}. (Default: 'phase')", default='phase')
    parser.add_option('--post', dest='postfit', action='store_true', help="Show postfit residuals. (Default: Don't show postfit.)", default=False)
    parser.add_option('--pre', dest='prefit', action='store_true', help="Show prefit residuals. (Default: Don't show prefit.)", default=False)
    parser.add_option('-l', '--legend', dest='legend', action='store_true', help="Show legend of frequencies. (Default: Do not show legend.)", default=False)
    parser.add_option('--mark-peri', dest='mark_peri', action='store_true', help="Mark passage of periastron. (Default: don't mark periastron.)", default=False)
    parser.add_option('--non-interactive', dest='interactive', action='store_false', help="Save figure and exit. (Default: Show plot, only save if requested.)", default=True)
    
    main()


