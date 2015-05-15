#!/usr/bin/env python

"""
pulse_energy_distribution.py

Calculate the energy of many Pulse objects and produce
a pulse energy distribution plot.
"""

import glob
import optparse
import sys
import warnings
import os.path

import numpy as np
import matplotlib.pyplot as plt

import pulse


def myhist(data, bins=50, *args, **kwargs):
    n, binedges = np.histogram(data, bins, new=True)
    binedges = binedges.repeat(2)
    n = np.concatenate(([0], n.repeat(2), [0]))
    n = np.clip(n, 0.1, n.max())
    
    plt.plot(binedges, n, *args, **kwargs)


def main():
    # Files can be provided as arguments, with "-g/--glob" option
    # or be contained in a file provided with "-f/--file" option.
    filenames = args + glob.glob(options.glob)
    if options.file is not None:
        if not os.path.exists(options.file):
            raise ValueError("File %s does not exist" % options.file)
        else:
            f = open(options.file, 'r')
            filenames += [fn.strip() for fn in f.readlines()]
    
    if not options.quiet:
        print "Number of files to consider: %d" % len(filenames)

    # Initialise lists to store energies
    on_energies = []
    off_energies = []

    # Collect energies
    for fn in filenames:
        if not os.path.exists(fn):
            continue
        prof = pulse.read_pulse_from_file(fn)
        (on, off) = prof.get_pulse_energies()
        on_energies.append(on)
        off_energies.append(off)
        
    on_mean = np.mean(on_energies)
    if not options.quiet:
        print "Average on-pulse energy: %f" % on_mean
    on = on_energies/on_mean
    off = off_energies/on_mean
    
    warnings.warn("Only plotting values with E/<E> > -5")
    on = on[on>-5]
    
    if not options.quiet:
        print "Number of pulses being plotted: %d" % len(on)
    
    # Plot results
    fig = plt.figure()
    myhist(on, bins=options.numbins, \
                color='k', linestyle='-', label="On Pulse")
    plt.xlabel("E/<E>")
    plt.ylabel("Number of Pulses")
    ymin, ymax = plt.ylim()
    plt.yscale("log")
    plt.ylim(0.5, ymax*2)
    plt.title(options.title)
    plt.legend(loc="best")
    if options.annotate:
        fig.text(0.05, 0.02, "Total # pulses plotted: %d" % on.size, \
                    ha='left', va='center', size='small')
    plt.savefig(options.savefn)
    if options.interactive:
        plt.show()


if __name__ == '__main__':
    parser = optparse.OptionParser(usage="%prog [options] pulse_files", \
                        description="Calculate the energy of many Pulse " \
                                    "objects and produce a pulse energy " \
                                    "distribution plot. Written by " \
                                    "Patrick Lazarus.", \
                        version="%prog v0.9 (by Patrick Lazarus)", \
                        prog="pulse_energy_distribution.py")
    parser.add_option('--debug', dest='debug', action='store_true', \
                        help="Display debugging information. (Default: " \
                                "don't display debugging information).", \
                        default=False)
    parser.add_option('-q', '--quiet', dest='quiet', action='store_true', \
                        help="Output less information. (Default: Output " \
                                "full information).", default=False)
    parser.add_option('-i', '--interactive', dest='interactive', \
                        action='store_true', default=False, \
                        help="Interactive mode. Show plot. (Default: Don't.)")
    parser.add_option('-a', '--annotate', dest='annotate', \
                        action='store_true', default=False, \
                        help="Add a note on the plot with some extra information. " \
                                "(Default: Don't.)")
    parser.add_option('-g', '--glob', dest='glob', help="Shell-style pattern " \
                        "defining files containing pulses. Be sure to " \
                        "use quotes.", default="")
    parser.add_option('-f', '--file', dest='file', help="File containing " \
                        "list of pulse file.", default=None)
    parser.add_option('-t', '--title', dest='title', help="Title on plot. " \
                        "(Default: No title.)", default="")
    parser.add_option('-s', '--savefn', dest='savefn', help="Filename to " \
                        "save plot as (Default: pulse_energy_distribution.ps)", \
                        default="pulse_energy_distribution.ps")
    parser.add_option('-n', '--numbins', dest='numbins', type='int', \
                        help="Number of bins in histogram. (Default: 50)", \
                        default=50)
    options, args = parser.parse_args()
    main()
