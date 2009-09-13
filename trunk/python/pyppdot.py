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

    def __str__(self):
        bfield, age, edot = self.get_computed_params()
        strings = ["PSR %s:" % self.name, \
                 "\tRA (J2000): %s, Dec (J2000): %s" % (self.raj, self.decj)]

        if self.p is not None:
            strings.append("\tPeriod (s): %f" % self.p)
        else:
            strings.append("\tPeriod (s): Not Measured")
            
        if self.pdot is not None:
            strings.append("\tP-dot (s/s): %0.3g" % self.pdot)
        else:
            strings.append("\tP-dot (s/s): Not Measured")
            
        if self.pdot is not None and self.pdot is not None:
            strings.extend(["\tB-field (G): %0.3g" % bfield, \
                                    "\tAge (yr): %0.3g" % age, \
                                    "\tE-dot (erg/s): %0.3g" % edot])
        return '\n'.join(strings)


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


def parse_pulsar_file(psrfn='pulsars.txt'):
    """Parse list of pulsars.
        Return a list of Pulsar objects.
    """
    pulsars = []
    psrfile = open(psrfn, 'r')
    for line in psrfile.readlines():
        split_line = line.split()
        name = split_line[1]
        raj = split_line[2]
        decj = split_line[3]
        if split_line[4] == '*':
            p = None
        else:
            p = float(split_line[4])
        if split_line[5] == '*':
            pdot = None
        else:
            pdot = float(split_line[5])
        if split_line[8] == '*':
            dm = None
        else:
            dm = float(split_line[8])
        if split_line[9] == '*':
            binary = False
        else:
            binary = True
        if split_line[10] == '*':
            assoc = False
        else:
            assoc = split_line[11]
        if split_line[11] == '*':
            psrtype = 'Radio'
        else:
            psrtype = split_line[11]
        pulsars.append(Pulsar(name, p, pdot, raj, decj, dm, \
                                binary, assoc, psrtype))
    return pulsars


def main():
    pulsars = parse_pulsar_file()
    for psr in pulsars[516:521]:
        print psr

if __name__=='__main__':
    main()
