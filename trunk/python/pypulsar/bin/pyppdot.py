#!/usr/bin/env python

"""
An interactive P-Pdot diagram plotter written in python
using matplotlib and numpy

          Patrick Lazarus, Sept 13th, 2009
"""

import optparse
import sys
import re
import os
import os.path
import types
import warnings

import matplotlib.pyplot as plt
import numpy as np

import psr_utils
from pypulsar.utils import colour

MARKER_OPTIONS = {'facecolor':'none', 'zorder':2, 'alpha':0.8, 'lw':2, 's':50} 
BINARY_MARKER = {'marker':'o', 'edgecolor':'g', 'label':'binary'}
RRAT_MARKER = {'marker':'s', 'edgecolor':'c', 'label':'rrat'}
MAGNETAR_MARKER = {'marker':'^', 'edgecolor':'m', 'label':'magnetar'}
SGR_MARKER = {'marker':'<', 'edgecolor':"#AB82FF", 'label':'sgr'}
AXP_MARKER = {'marker':'>', 'edgecolor':"#FF7F24", 'label':'axp'}
SNR_MARKER = {'marker':(4,1,0), 'edgecolor':'y', 'label':'snr'}

# Picker determines how far mouse can be from point to select it
PICKER = 100

class Pulsar:
    def __init__(self, name, p, pdot, raj, decj, dm, \
                    binarytype, assoc, psrtype):
        """Pulsar object.
            Attributes:
                name - Name
                p - Period
                pdot - Period derivative
                raj - RA J2000
                decj - DEC J2000
                dm - DM
                binarytype - Binary type
                assoc - Associations
                psrtype - Pulsar type
        """
        self.name = name
        self.p = p
        self.pdot = pdot
        self.raj = raj
        self.decj = decj
        self.dm = dm
        self.binarytype = binarytype
        self.assoc = assoc
        self.psrtype = psrtype
        self.rrat = (psrtype is not None and psrtype!='No info' and 'rrat' in psrtype.lower())
        self.magnetar = (psrtype is not None and psrtype!='No info' and \
                            ('axp' in psrtype.lower() or 'sgr' in psrtype.lower()))
        self.sgr = (psrtype is not None and psrtype!='No info' and \
                            ('axp' in psrtype.lower() and 'sgr' in assoc.lower()))
        self.axp = (psrtype is not None and psrtype!='No info' and \
                            ('axp' in psrtype.lower() and 'sgr' not in assoc.lower()))
        self.snr = (assoc is not None and assoc!='No info' and 'snr' in assoc.lower())
        self.binary = (binarytype is not None and binarytype!='No info')

    def get_computed_params(self):
        """Return pulsar's B-field, age and Edot in a tuple.
            If either P or Pdot is not known, return tuple of Nones
        """
        return params_from_ppdot(self.p, self.pdot)

    def get_info(self, extended=False):
        """Return information string describing pulsar.
            If extended is True then return all available information.
        """
        bfield, age, edot = self.get_computed_params()
        strings = [colour.cstring("PSR %s" % self.name, \
                                    fg='default', bg='default', underline=1, bold=1), \
                 "\tRA (J2000): %s, Dec (J2000): %s" % (self.raj, self.decj)]

        if self.p is not None:
            strings.append("\tPeriod (s): %f" % self.p)
        else:
            strings.append("\tPeriod (s): Not Measured")
            
        if self.pdot is not None:
            strings[-1] += ", P-dot (s/s): %0.3g" % self.pdot
        else:
            strings[-1] += ", P-dot (s/s): Not Measured"
            
        if self.pdot is not None and self.pdot is not None:
            strings.extend(["\tB-field (G): %0.3g" % bfield, \
                                    "\tAge (%s): %0.3g" % units_age(age), \
                                    "\tE-dot (erg/s): %0.3g" % edot])
        if extended:
            strings.extend(["\tBinary type: %s" % self.binarytype, \
                            "\tAssociations: %s" % self.assoc, \
                            "\tPulsar type: %s" % self.psrtype])
        return '\n'.join(strings)

    
    def __str__(self):
        return self.get_info(extended=False)


def units_age(age):
    prefix = ["", "k", "M", "G"]
    m = int(np.log10(age)/3)
    if m >= len(prefix):
        m = len(prefix)-1
    r = age/10**(m*3)
    return ("%syr" % prefix[m], r)


def pdot_from_edot(p, edot):
    """Return pdot (in s/s) that a pulsar with period, p (in s), would require 
        to have a particular observed edot (in erg/s).

        (Use eq. 3.6 from Lorimer and Kramer.)
    """
    # The multiplicative constant is: 1e-15/3.95e31
    return 2.5316455696202532e-47*edot*(p**3)


def p_from_edot(pdot, edot):
    """Return period (in s) that a pulsar with pdot (in s/s), would require 
        to have a particular observed edot (in erg/s).

        (Use eq. 3.6 from Lorimer and Kramer.)
    """
    # The multiplicative constant is: 1e-15/3.95e31
    return (pdot/(2.5316455696202532e-47*edot))**(1/3.0)


def pdot_from_bfield(p, bfield):
    """Return pdot (in s/s) that a pulsar with period, p (in s), would require
        to have a particular observed B-field (in G).

        (Use eq. 3.15 from Lorimer and Kramer.)
    """
    # The multiplicative constant is: 1e-15/(1e12)**2
    return 1e-39*(bfield)**2/p


def p_from_bfield(pdot, bfield):
    """Return period (in s) that a pulsar with pdot (in s/s), would require
        to have a particular observed B-field (in G).

        (Use eq. 3.15 from Lorimer and Kramer.)
    """
    # The multiplicative constant is: 1e-15/(1e12)**2
    return 1e-39*(bfield)**2/pdot


def pdot_from_age(p, age):
    """Return pdot (in s/s) that a pulsar with period, p (in s), would require
        to have a particular observed age (in yrs).

        (Use eq. 3.12 from Lorimer and Kramer.)
    """
    return p/age/(2.0*psr_utils.SECPERJULYR)


def p_from_age(pdot, age):
    """Return period (in s) that a pulsar with pdot (in s/s), would require
        to have a particular observed age (in yrs).

        (Use eq. 3.12 from Lorimer and Kramer.)
    """
    return pdot*age*(2.0*psr_utils.SECPERJULYR)


def params_from_ppdot(p, pdot):
    """Return B-field, age and Edot in a tuple given p and pdot.
        If either P or Pdot is None, return tuple of Nones
    """
    if p is None or pdot is None:
        bfield = None
        age = None
        edot = None
    else:
        f, fdot = psr_utils.p_to_f(p, pdot)
        bfield = psr_utils.pulsar_B(f, fdot)
        age = psr_utils.pulsar_age(f, fdot)
        edot = psr_utils.pulsar_edot(f, fdot)
    return (bfield, age, edot)


def plot_data(pulsars, hightlight=[], binaries=False, rrats=False, \
                magnetars=False, snrs=False, axp=False, sgr=False, \
                edots=[], ages=[], bsurfs=[]):
    """Plot P-Pdot diagram using list of pulsars provided.
        binaries - boolean, initially plot binary markers
        rrats - boolean, initially plot RRAT markers
        magnetars - boolean, initially plot magnetar markers
        snrs - boolean, initially plot SNR markers
        edots - values for lines of constant edot
        ages - values for lines of constant ages
        bsurfs - values for lines of constant surfave b-field
    """
    global scatt_psrs, scatt_psrs_hl
    periods = np.array([x.p for x in pulsars \
                        if x.p is not None and x.pdot is not None])
    pdots = np.array([x.pdot for x in pulsars \
                        if x.p is not None and x.pdot is not None])
    
    ax = plt.axes()
    scatt_psrs = ax.scatter(periods, pdots, c='k', s=6, \
                            label='pulsars', picker=PICKER, zorder=0)

    # Pulsars to highlight
    numhl = 0
    if len(highlight):
        periods_hl = np.array([h.p for h in highlight \
                            if h.p is not None and h.pdot is not None])
        pdots_hl = np.array([h.pdot for h in highlight \
                            if h.p is not None and h.pdot is not None])
        scatt_psrs_hl = ax.scatter(periods_hl, pdots_hl, \
                                    c='r', s=50, label='highlight', picker=PICKER, \
                                    marker=(5,1,0), edgecolors='r', zorder=1)
        numhl = periods_hl.size

    # Mark binaries
    periods_bnry = np.array([x.p for x in pulsars \
                            if x.p is not None and x.pdot is not None \
                                and x.binary==True])
    pdots_bnry = np.array([x.pdot for x in pulsars \
                            if x.p is not None and x.pdot is not None \
                                and x.binary==True])
    global scatt_binaries
    if periods_bnry.size:
        scatter_options = MARKER_OPTIONS.copy()
        scatter_options.update(BINARY_MARKER)
        scatt_binaries = ax.scatter(periods_bnry, pdots_bnry, \
                                    **scatter_options)
        if not binaries:
            # Hide binaries for now
            scatt_binaries.set_visible(False)
    else:
        scatt_binaries = None

    # Mark RRATs
    periods_rrat = np.array([x.p for x in pulsars+highlight \
                            if x.p is not None and x.pdot is not None \
                                and x.rrat==True])
    pdots_rrat = np.array([x.pdot for x in pulsars+highlight \
                            if x.p is not None and x.pdot is not None \
                                and x.rrat==True])
    global scatt_rrats
    if periods_rrat.size:
        scatter_options = MARKER_OPTIONS.copy()
        scatter_options.update(RRAT_MARKER)
        scatt_rrats = ax.scatter(periods_rrat, pdots_rrat, \
                                    **scatter_options)
        if not rrats:
            # Hide RRATs for now
            scatt_rrats.set_visible(False)
    else:
        scatt_rrats = None

    # Mark magnetars
    periods_mag = np.array([x.p for x in pulsars \
                            if x.p is not None and x.pdot is not None \
                                and x.magnetar==True])
    pdots_mag = np.array([x.pdot for x in pulsars \
                            if x.p is not None and x.pdot is not None \
                                and x.magnetar==True])
    global scatt_magnetars
    if periods_mag.size:
        scatter_options = MARKER_OPTIONS.copy()
        scatter_options.update(MAGNETAR_MARKER)
        scatt_magnetars = ax.scatter(periods_mag, pdots_mag, \
                                    **scatter_options)
        if not magnetars:
            # Hide magnetars for now
            scatt_magnetars.set_visible(False)
    else:
        scatt_magnetars = None

    # Mark sgrs
    periods_sgr = np.array([x.p for x in pulsars \
                            if x.p is not None and x.pdot is not None \
                                and x.sgr==True])
    pdots_sgr = np.array([x.pdot for x in pulsars \
                            if x.p is not None and x.pdot is not None \
                                and x.sgr==True])
    global scatt_sgr
    if periods_sgr.size:
        scatter_options = MARKER_OPTIONS.copy()
        scatter_options.update(SGR_MARKER)
        scatt_sgr = ax.scatter(periods_sgr, pdots_sgr, \
                                    **scatter_options)
        if not sgr:
            # Hide sgr for now
            scatt_sgr.set_visible(False)
    else:
        scatt_sgr = None

    # Mark axps
    periods_axp = np.array([x.p for x in pulsars \
                            if x.p is not None and x.pdot is not None \
                                and x.axp==True])
    pdots_axp = np.array([x.pdot for x in pulsars \
                            if x.p is not None and x.pdot is not None \
                                and x.axp==True])
    global scatt_axp
    if periods_axp.size:
        scatter_options = MARKER_OPTIONS.copy()
        scatter_options.update(AXP_MARKER)
        scatt_axp = ax.scatter(periods_axp, pdots_axp, \
                                    **scatter_options)
        if not axp:
            # Hide axp for now
            scatt_axp.set_visible(False)
    else:
        scatt_axp = None

    # Mark SNRs
    periods_snr = np.array([x.p for x in pulsars \
                            if x.p is not None and x.pdot is not None \
                                and x.snr==True])
    pdots_snr = np.array([x.pdot for x in pulsars \
                            if x.p is not None and x.pdot is not None \
                                and x.snr==True])
    global scatt_snrs
    if periods_snr.size:
        scatter_options = MARKER_OPTIONS.copy()
        scatter_options.update(SNR_MARKER)
        scatt_snrs = ax.scatter(periods_snr, pdots_snr, \
                                    **scatter_options)
        if not snrs:
            # Hide SNRs for now
            scatt_snrs.set_visible(False)
    else:
        scatt_snrs = None

    plt.xlabel("Period (s)")
    plt.ylabel(r"$\mathsf{\dot P}$ (s/s)")
    plt.title(r"$\mathsf{P-\dot P}$ Diagram")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim((0.001, 100))
    plt.ylim((10**-22, 10**-8))
    draw_lines(bsurfs, edots, ages)

    print "Plot Inventory:"
    print "\tNumber of pulsars:", (periods.size + numhl)
    print "\tNumber highlighted:", numhl
    print "\tNumber of RRATs:", periods_rrat.size
    print "\tNumber of magnetars:", periods_mag.size
    print "\tNumber of binaries:", periods_bnry.size
    print "\tNumber of SNR associations:", periods_snr.size

def draw_lines(bfields=[], edots=[], ages=[], label=False):
    """Draw multiple lines of constant B-field, edot and age.
    """
    plimits = plt.xlim()
    pdotlimits = plt.ylim()
    periods = np.logspace(np.log10(plimits[0]), np.log10(plimits[1]), 100)
    pdots = np.logspace(np.log10(pdotlimits[0]), np.log10(pdotlimits[1]), 100)
    bfield_vals = np.zeros((periods.size, pdots.size))
    edot_vals = np.zeros((periods.size, pdots.size))
    age_vals = np.zeros((periods.size, pdots.size))
    for i in np.arange(periods.size):
        for j in np.arange(pdots.size):
            f, fdot = psr_utils.p_to_f(periods[i], pdots[j])
            # For some reason the transpose is expected by 'contour'...
            bfield_vals[j][i] = psr_utils.pulsar_B(f, fdot)
            edot_vals[j][i] = psr_utils.pulsar_edot(f, fdot)
            age_vals[j][i] = psr_utils.pulsar_age(f, fdot)
    if len(bfields) > 1:
        cs_bfield = plt.contour(periods, pdots, bfield_vals, bfields, \
                                colors='k', linestyles='dashdot')
        if label: plt.clabel(cs_bfield, fmt=r'%g $G$')
    if len(edots) > 1:
        cs_edot = plt.contour(periods, pdots, edot_vals, edots, \
                                colors='k', linestyles='dashed')
        if label: plt.clabel(cs_edot, fmt=r'%g $erg/s$')
    if len(ages) > 1:
        cs_age = plt.contour(periods, pdots, age_vals, ages, \
                                colors='k', linestyles='dotted')
        if label: plt.clabel(cs_age, fmt=r'%g $yr$')
    plt.xlim(plimits)
    plt.ylim(pdotlimits)
    

def quit():
    print "Quiting..."
    sys.exit(0)


def savefigure(savefn='./ppdot.png'):
    print "Saving plot to %s" % savefn
    plt.savefig(savefn, orientation='landscape', papertype='letter')


def mousepress(event):
    """Event handler for MouseEvent ('button_press_event').
    """
    if event.inaxes and event.button == 2:
        p = event.xdata
        pdot = event.ydata
        bfield, age, edot = params_from_ppdot(p, pdot)
        print "Coordinates:"
        print "\tPeriod (s): %g, P-dot (s/s): %g" % (p, pdot)
        print "\tB-field (G): %g" % bfield
        print "\tAge (yr): %g" % age
        print "\tE-dot (erg/s): %g" % edot


def pick(event):
    """Event handler for PickEvent.
        Display pulsar information.
    """
    global pulsars
    global highlight
    if event.mouseevent.button == 1 and len(event.ind) == 1:
        index = event.ind[0]
        sys.stdout.write("Pulsar selected: ")
        if event.artist.get_label() == 'pulsars':
            selected = pulsars[index]
        elif event.artist.get_label() == 'highlight':
            selected = highlight[index]
        else:
            print "What was selected?! Error?"
            return
        if event.mouseevent.key == 'shift':
            # provide complete information
            print selected.get_info(extended=True)
        else:
            print selected.get_info(extended=False)


def keypress(event):
    """Event handler for KeyEvent.
    """
    if type(event.key) == types.StringType:
        
        if event.key.lower() == 'q':
            # Quit program
            quit()
        elif event.key.lower() == 's':
            # Save figure
            savefigure()
        elif event.key.lower() == 'z':
            # Toggle zoom mode
            print "Toggling zoom mode..."
            event.canvas.toolbar.zoom()
        elif event.key.lower() == 'o':
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
        elif event.key.lower() == 'b':
            # Mark binaries
            global scatt_binaries
            if scatt_binaries is not None:
                print "Toggling binaries..."
                visible = scatt_binaries.get_visible()
                # visible is True/False. 'not visible' will toggle state.
                scatt_binaries.set_visible(not visible)
                event.canvas.draw()
            else:
                print "(No binaries)"
        elif event.key.lower() == 'r':
            # Mark RRATs
            global scatt_rrats
            if scatt_rrats is not None:
                print "Toggling RRATs..."
                visible = scatt_rrats.get_visible()
                # visible is True/False. 'not visible' will toggle state.
                scatt_rrats.set_visible(not visible)
                event.canvas.draw()
            else:
                print "(No RRATs)"
        elif event.key.lower() == 'm':
            # Mark magnetars
            global scatt_magnetars
            if scatt_magnetars is not None:
                print "Toggling magnetars..."
                visible = scatt_magnetars.get_visible()
                # visible is True/False. 'not visible' will toggle state.
                scatt_magnetars.set_visible(not visible)
                event.canvas.draw()
            else:
                print "(No magnetars)"
        elif event.key.lower() == 'g':
            # Mark sgrs
            global scatt_sgr
            if scatt_sgr is not None:
                print "Toggling SGRs..."
                visible = scatt_sgr.get_visible()
                # visible is True/False. 'not visible' will toggle state.
                scatt_sgr.set_visible(not visible)
                event.canvas.draw()
            else:
                print "(No SGRs)"
        elif event.key.lower() == 'a':
            # Mark axps
            global scatt_axp
            if scatt_axp is not None:
                print "Toggling AXPs..."
                visible = scatt_axp.get_visible()
                # visible is True/False. 'not visible' will toggle state.
                scatt_axp.set_visible(not visible)
                event.canvas.draw()
            else:
                print "(No AXPs)"
        elif event.key.lower() == 'n':
            # Mark SNRs
            global scatt_snrs
            if scatt_snrs is not None:
                print "Toggling SNR associations..."
                visible = scatt_snrs.get_visible()
                # visible is True/False. 'not visible' will toggle state.
                scatt_snrs.set_visible(not visible)
                event.canvas.draw()
            else:
                print "(No SNR associations)"
        elif event.key == 'h':
            # Display help
            print "Helping..."
            print "-"*80
            print "Help - Hotkeys definitions:"
            print "\th - Display this help"
            print "\tq - Quit"
            print "\ts - Save current plot to PostScript file"
            print "\tz - Toggle Zoom-mode on/off"
            print "\to - Go to original view"
            print "\t< - Go to previous view"
            print "\t> - Go to next view"
            print "\tb - Toggle binary marking"
            print "\tr - Toggle RRAT marking"
            print "\tm - Toggle magnetar marking"
            print "\tn - Toggle SNR marking"
            print "\t[Left mouse] - Select pulsar (display info in terminal)"
            print "\t             - Select zoom region (if Zoom-mode is on)"
            print "\t[Middle mouse] - Display P, P-dot, B, E-dot and age at mouse pointer"
            print "-"*80
            

def create_plot(pulsars, highlight=[], interactive=True, **kwargs):
    """Create the plot and set up event handlers.
        pulsars - a list of pulsar objects
        highlight - a list of pulsar objects to highlight
        interactive - boolean, enter interative mode after making the plot
        binaries - boolean, intial state for marking binaries on the plot
        rrats - boolean, intial state for marking rrats on the plot
        magnetars - boolean, intial state for marking magnetars on the plot
        snrs - boolean, intial state for marking snrs on the plot
    """
    fig = plt.figure(figsize=(11,8.5))
    fig.canvas.set_window_title("P-Pdot")
    plot_data(pulsars, highlight, **kwargs)
    
    if interactive:
        # Before setting up our own event handlers delete matplotlib's
        # default 'key_press_event' handler.
        defcids = fig.canvas.callbacks.callbacks['key_press_event'].keys()
        for cid in defcids:
            fig.canvas.callbacks.disconnect(cid)
        
        # Register event callbacks function and show the plot
        cid_keypress = fig.canvas.mpl_connect('key_press_event', keypress)
        cid_mousepress = fig.canvas.mpl_connect('button_press_event', mousepress)
        cid_pick = fig.canvas.mpl_connect('pick_event', pick)
        plt.ion()
        plt.show()


def parse_pulsar_file(psrfn='pulsars.txt', indent=""):
    """Parse list of pulsars.
        Return a list of Pulsar objects.
    """
    print indent+"Parsing file (%s)" % psrfn
    pulsars = []
    nonplottable = 0
    
    # Check if file exists
    if not os.path.exists(psrfn):
        print indent+"    File not found: %s" % psrfn
        return pulsars
    
    psrfile = open(psrfn, 'r')
    for line in psrfile.readlines():
        line = line.partition('#')[0].strip()
        if len(line) == 0:
            continue
        split_line = line.split()

        if split_line[0].upper() == 'INCLUDE':
            # include pulsars from other files.
            for fn in split_line[1:]:
                print indent+"    INCLUDE'ing another file."
                newindent=indent+'    '
                dir = os.path.split(psrfn)[0]
                pulsars += parse_pulsar_file(os.path.join(dir,fn), indent=newindent)
            continue
        else:
            name = split_line[0]
        
        if split_line[1] == '*':
            p = None
            # No period measured, can't be plotted
            nonplottable += 1
            continue
        else:
            p = float(split_line[1])
            
        if split_line[2] == '*':
            pdot = None
            # No pdot measured, can't be plotted
            nonplottable += 1
            continue
        else:
            pdot = float(split_line[2])
            
        if len(split_line) < 4 or split_line[3] == '*':
            raj = None
        else:    
            raj = split_line[3]
            
        if len(split_line) < 5 or split_line[4] == '*':
            decj = None
        else:
            decj = split_line[4]
            
        if len(split_line) < 6 or split_line[5] == '*':
            dm = None
        else:
            dm = float(split_line[5])
            
        if len(split_line) < 7:
            binarytype = "No info"
        elif split_line[6] == '*':
            binarytype = None
        else:
            binarytype = split_line[6]
            
        if len(split_line) < 8:
            assoc = "No info"
        elif split_line[7] == '*':
            assoc = None
        else:
            assoc = split_line[7]
            
        if len(split_line) < 9:
            psrtype = "No info"
        elif split_line[8] == '*':
            psrtype = 'Radio'
        else:
            psrtype = split_line[8]
        pulsars.append(Pulsar(name, p, pdot, raj, decj, dm, \
                                binarytype, assoc, psrtype))
    print indent+"    Number of pulsars that cannot be plotted (no P or Pdot): %d" % nonplottable
    print indent+"    Done parsing file."
    return pulsars


def parse_options():
    (options, sys.argv) = parser.parse_args()
    if options.def_lines:
        options.edots = [1e30, 1e33, 1e36]
        options.bsurfs = [1e10, 1e12, 1e14]
        options.ages = [1e3, 1e6, 1e9]
    return options


def main():
    global pulsars
    global highlight
    options = parse_options()
    highlight = []
    if options.files:
        pulsars = []
        for file in options.files:
            pulsars += parse_pulsar_file(file)
    else:
        pulsars = parse_pulsar_file()
    for hl in options.highlight:
        highlight += parse_pulsar_file(hl)

    # Eliminate duplicate pulsars, based on name
    psr_dict = {}
    for psr in pulsars:
        psr_dict[psr.name] = psr
    for hl in highlight:
        if hl.name in psr_dict:
            del psr_dict[hl.name]
    pulsars = psr_dict.values()
    if len(pulsars)+len(highlight):
        create_plot(pulsars, highlight, binaries=options.binaries, \
                magnetars=options.magnetars, rrats=options.rrats, \
                snrs=options.snrs, edots=options.edots, ages=options.ages, \
                bsurfs=options.bsurfs)


if __name__=='__main__':
    parser = optparse.OptionParser()
    parser.add_option('-f', '--file', dest='files', type='string', \
                        action='append', help="File containing a list "
                        "of pulsars to display with ATNF catalogue. "
                        "Each pulsar should be on a separate row with "
                        "the following format:\nName period pdot dm "
                        "binary associations pulsar_type."
                        "\nEach column should contain a single string "
                        "(no space), and '*' should be used as a null value.", \
                        default=[])
    parser.add_option('--highlight', dest='highlight', type='string', \
                        action='append', help="File containing a list "
                        "of pulsars to display with ATNF catalogue. "
                        "These pulsars will be highlighed (displayed "
                        "with a star instead of a point). See -f/--file "
                        "option for formatting.", default=[])
    parser.add_option('-e', '--edot', dest='edots', type='float', \
                        action='append', help="Value, in erg/s, to plot "
                        "a line of constant E-dot. Multiple -e/--edot options "
                        "can be provided.", default=[])
    parser.add_option('-a', '--age', dest='ages', type='float', \
                        action='append', help="Value, in yr, to plot "
                        "a line of constant age. Multiple -a/--age options "
                        "can be provided.", default=[])
    parser.add_option('-b', '--bsurf', dest='bsurfs', type='float', \
                        action='append', help="Value, in G, to plot "
                        "a line of constant surface B-field. Multiple "
                        "-b/--bsurf options can be provided.", \
                        default=[])
    parser.add_option('--def-lines', dest='def_lines', action='store_true', \
                        help="Plot default lines\n"
                        "E-dot (erg/s): 1e30, 1e33, 1e36\n"
                        "B-field (G): 1e10, 1e12, 1e14\n"
                        "Age (yr): 1e3, 1e6, 1e9", default=False)
    parser.add_option('--binaries', dest='binaries', action='store_true', \
                        help="Mark binary pulsars. This is the initial state, "
                        "binary marking can be toggled interactively. "
                        "(Default: Don't distinguish binaries.)", \
                        default=False)
    parser.add_option('--rrats', dest='rrats', action='store_true', \
                        help="Mark RRATs. This is the initial state, RRAT "
                        "marking can be toggled interactively. (Default: "
                        "Don't distinguish RRATs.)", default=False)
    parser.add_option('--magnetars', dest='magnetars', \
                        action='store_true', help="Mark magnetars. This "
                        "is the initial state, magnetar marking can be "
                        "toggled interactively. (Default: Don't distinguish "
                        "magnetars.)", default=False)
    parser.add_option('--snrs', dest='snrs', action='store_true', \
                        help="Mark supernova remnant associations. This "
                        "is the initial state, SNR marking can be toggled "
                        "interactively. (Default: Don't distinguish SNR "
                        "associations.)", default=False)
    main()
