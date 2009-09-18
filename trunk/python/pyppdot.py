#!/usr/bin/env python

# An interactive P-Pdot diagram plotter written in python
# using matplotlib and numpy
#
#           Patrick Lazarus, Sept 13th, 2009

import optparse
import sys
import re
import os
import types
import warnings

import matplotlib.pyplot as plt
import numpy as np

import psr_utils
import colour

class Pulsar:
    def __init__(self, name, p, pdot, raj, decj, dm, binary, assoc, psrtype):
        """Pulsar object.
            Attributes:
                name - Name
                p - Period
                pdot - Period derivative
                raj - RA J2000
                decj - DEC J2000
                dm - DM
                binary - Binary
                assoc - Associations
                psrtype - Pulsar type
        """
        self.name = name
        self.p = p
        self.pdot = pdot
        self.raj = raj
        self.decj = decj
        self.dm = dm
        self.binary = binary
        self.assoc = assoc
        self.psrtype = psrtype

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
                                    "\tAge (yr): %0.3g" % age, \
                                    "\tE-dot (erg/s): %0.3g" % edot])
        if extended:
            strings.extend(["\tBinary: %s" % self.binary, \
                            "\tAssociations: %s" % self.assoc, \
                            "\tPulsar type: %s" % self.psrtype])
        return '\n'.join(strings)

    def __str__(self):
        return self.get_info(extended=False)


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


def plot_data(pulsars, hightlight=[]):
    """Plot P-Pdot diagram using list of pulsars provided.
    """
    global scatt_psrs, scatt_psrs_hl
    periods = np.array([x.p for x in pulsars \
                        if x.p is not None and x.pdot is not None])
    pdots = np.array([x.pdot for x in pulsars \
                        if x.p is not None and x.pdot is not None])
    ax = plt.axes()
    scatt_psrs = ax.scatter(np.log10(periods), np.log10(pdots), c='k', s=6, \
                            label='pulsars', picker=10)

    # Pulsars to highlight
    if len(highlight):
        print "Highlighting %d pulsars." % len(highlight)
        periods_hl = np.array([h.p for h in highlight \
                            if h.p is not None and h.pdot is not None])
        pdots_hl = np.array([h.pdot for h in highlight \
                            if h.p is not None and h.pdot is not None])
        scatt_psrs_hl = ax.scatter(np.log10(periods_hl), np.log10(pdots_hl), \
                                    c='r', s=50, label='highlight', picker=10, \
                                    marker=(5,1,0), edgecolors='none')

    # Mark binaries
    periods_bnry = np.array([x.p for x in pulsars \
                            if x.p is not None and h.pdot is not None \
                                and x.binary==True])
    pdots_bnry = np.array([x.pdot for x in pulsars \
                            if x.p is not None and h.pdot is not None \
                                and x.binary==True])
    print "%d binary pulsars." % periods_bnry.size
    global scatt_binaries
    scatt_binaries = ax.scatter(np.log10(periods_bnry), np.log10(pdots_bnry), \
                                s=50, lw=2, facecolor='none', \
                                edgecolor='g', label='binary')
    if True:
        # Hide binaries for now
        scatt_binaries.set_visible(False)
    plt.xlabel("log Period")
    plt.ylabel("log P-dot")
    plt.title("P vs. P-dot")
    plimits = np.array((0.001, 100))
    pdotlimits = np.array((10**-22, 10**-8))
    draw_bfield_line(1e10, plimits, pdotlimits)
    draw_bfield_line(1e12, plimits, pdotlimits)
    draw_bfield_line(1e14, plimits, pdotlimits)
    draw_edot_line(1e30, plimits, pdotlimits)
    draw_edot_line(1e33, plimits, pdotlimits)
    draw_edot_line(1e36, plimits, pdotlimits)
    draw_age_line(1e3, plimits, pdotlimits)
    draw_age_line(1e6, plimits, pdotlimits)
    draw_age_line(1e9, plimits, pdotlimits)
    

def draw_bfield_line(bfield, plimits, pdotlimits):
    """Draw and annotate a line of constant B-field.
    """
    ax = plt.gca()
    ax.plot(np.log10(plimits), np.log10(pdot_from_bfield(plimits, bfield)), 'k-.')
    ax.set_xlim(np.log10(plimits))
    ax.set_ylim(np.log10(pdotlimits))
    annotate_line(bfield, 'G', plimits, pdotlimits, pdot_from_bfield, p_from_bfield)


def draw_edot_line(edot, plimits, pdotlimits):
    """Draw and annotate a line of constant E-dot.
    """
    ax = plt.gca()
    ax.plot(np.log10(plimits), np.log10(pdot_from_edot(plimits, edot)), 'k--')
    ax.set_xlim(np.log10(plimits))
    ax.set_ylim(np.log10(pdotlimits))
    annotate_line(edot, 'erg/s', plimits, pdotlimits, pdot_from_edot, p_from_edot)


def draw_age_line(age, plimits, pdotlimits):
    """Draw and annotate a line of constant age.
    """
    ax = plt.gca()
    ax.plot(np.log10(plimits), np.log10(pdot_from_age(plimits, age)), 'k:')
    ax.set_xlim(np.log10(plimits))
    ax.set_ylim(np.log10(pdotlimits))
    annotate_line(age, 'yr', plimits, pdotlimits, pdot_from_age, p_from_age)


def annotate_line(value, units, xlimits, ylimits, pdot_from_value, p_from_value):
    """Neatly annotate a line of constant B-field/E-dot/age on figure.
        value: the value of the line
        units: the units to display
        xlimits, ylimits: min/max values of x- and y- axes
        pdot_from_value: a function to calculate pdot given p and value
        p_from_value: a function to calculate p given pdot and value
    """
    # Generate text for annotation
    exponent = int(np.log10(value))
    # Round matissa to 1 decimal place
    mantissa = np.round(10**(np.log10(value)-exponent),1)
    if mantissa != 1.0:
        mult = r"%0.1f \times " % mantissa
    else:
        mult = ""
    label_text = r"$%s10^{%d} %s$" % (mult, exponent, units)
    
    trans = plt.gca().transData + plt.gcf().transFigure.inverted()
    xi,xf = np.log10(xlimits)
    yi = np.log10(pdot_from_value(10**xi, value))
    yf = np.log10(pdot_from_value(10**xf, value))
    xi_trans,yi_trans = trans.transform((xi,yi))
    xf_trans,yf_trans = trans.transform((xf,yf))
    slope = float(yf_trans-yi_trans)/float(xf_trans-xi_trans)
    angle = np.arctan(slope)*180/np.pi
    ylim_i, ylim_f = plt.gca().get_ylim()
    if yf > ylim_i and yf < ylim_f:
        # line intersects right edge of plot window
        pass
    elif yi > ylim_i and yf < ylim_i:
        # line intersects bottom edge of plot window
        yf = np.log10(ylimits[0])
        xf = np.log10(p_from_value(10**yf, value))
        xf_trans,yf_trans = trans.transform((xf,yf))
    elif yi < ylim_f and yf > ylim_f:
        # line intersects top edge of plot window
        yf = np.log10(ylimits[1])
        xf = np.log10(p_from_value(10**yf, value))
        xf_trans,yf_trans = trans.transform((xf,yf))
    else:
        # line doesn't pass through plot window
        return
    # Write annotation on figure
    plt.figtext(xf_trans, yf_trans, label_text, \
                    rotation=angle, verticalalignment='center')
    

def quit():
    print "Quiting..."
    sys.exit(0)


def savefigure(savefn='./ppdot.ps'):
        print "Saving plot to %s" % savefn
        plt.savefig(savefn, orientation='landscape', papertype='letter')


def mousepress(event):
    """Event handler for MouseEvent ('button_press_event').
    """
    if event.inaxes and event.button == 2:
        p = 10**(event.xdata)
        pdot = 10**(event.ydata)
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
            print "Toggling binaries..."
            global scatt_binaries
            visible = scatt_binaries.get_visible()
            # visible is True/False. 'not visible' will toggle state.
            scatt_binaries.set_visible(not visible)
            event.canvas.draw()
        elif event.key == 'h':
            # Display help
            print "Helping..."
            print "-"*80
            print "Help - Hotkeys definitions:"
            print "\th - Display this help"
            print "\tq - Quit"
            print "\ts - Save current plot(s) to PostScript file"
            print "\tz - Toggle Zoom-mode on/off"
            print "\to - Go to original view"
            print "\t< - Go to previous view"
            print "\t> - Go to next view"
            print "\tb - Toggle binary marking"
            print "\t[Left mouse] - Select pulsar (display info in terminal)"
            print "\t             - Select zoom region (if Zoom-mode is on)"
            print "\t[Middle mouse] - Display P, P-dot, B, E-dot and age at mouse pointer"
            print "-"*80
            

def create_plot(pulsars, highlight=[], interactive=True):
    """Create the plot and set up event handlers.
    """
    fig = plt.figure(figsize=(11,8.5))
    plot_data(pulsars, highlight)
    
    # Register event callbacks function and show the plot
    cid_keypress = fig.canvas.mpl_connect('key_press_event', keypress)
    cid_mousepress = fig.canvas.mpl_connect('button_press_event', mousepress)
    cid_pick = fig.canvas.mpl_connect('pick_event', pick)
    
    if interactive:
        plt.ion()
        plt.show()


def parse_pulsar_file(psrfn='pulsars.txt'):
    """Parse list of pulsars.
        Return a list of Pulsar objects.
    """
    print "Parsing file (%s)" % psrfn
    pulsars = []
    nonplottable = 0
    psrfile = open(psrfn, 'r')
    for line in psrfile.readlines():
        split_line = line.split()
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
            binary = "No info"
        elif split_line[6] == '*':
            binary = False
        else:
            binary = True
            
        if len(split_line) < 8:
            assoc = "No info"
        elif split_line[7] == '*':
            assoc = False
        else:
            assoc = split_line[7]
            
        if len(split_line) < 9:
            psrtype = "No info"
        elif split_line[8] == '*':
            psrtype = 'Radio'
        else:
            psrtype = split_line[8]
        pulsars.append(Pulsar(name, p, pdot, raj, decj, dm, \
                                binary, assoc, psrtype))
    print "\tNumber of pulsars that cannot be plotted (no P or Pdot): %d" % nonplottable
    return pulsars


def parse_options():
    (options, sys.argv) = parser.parse_args()
    return options

def main():
    global pulsars
    global highlight
    options = parse_options()
    pulsars = parse_pulsar_file()
    highlight = []
    for file in options.files:
        pulsars += parse_pulsar_file(file)
    for hl in options.highlight:
        highlight += parse_pulsar_file(hl)
    create_plot(pulsars)

if __name__=='__main__':
    parser = optparse.OptionParser()
    parser.add_option('-f', '--file', dest='files', type='string', action='append', help="File containing a list of pulsars to display with ATNF catalogue. Each pulsar should be on a separate row with the following format:\nName period pdot dm binary associations pulsar_type.\nEach column should contain a single string (no space), and '*' should be used as a null value.", default=[])
    parser.add_option('--highlight', dest='highlight', type='string', action='append', help="File containing a list of pulsars to display with ATNF catalogue. These pulsars will be highlighed (displayed with a star instead of a point). See -f/--file option for formatting.", default=[])
    main()
