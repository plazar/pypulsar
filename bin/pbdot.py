#!/usr/bin/env python

"""
pbdot.py

A script to determine when Pb-dot should be detectable given binary system paramters
and the current uncertainty on Pb.

Uses equation 8.52 from Lorimer and Kramer's Pulsar Handbook.


Patrick Lazarus
July 5, 2012
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker

import psr_utils

TSUN = 4.925490947e-6 # Solar mass in units of time (see LK pg 217)

MP_MIN = 1.2 # Minimum pulsar mass (in Msun)
MP_MAX = 3.0 # Maximum pulsar mass (in Msun)
MC_MIN = 0.9 # Minimum pulsar mass (in Msun)
MC_MAX = 3.0 # Maximum pulsar mass (in Msun)

ORBITAL_PERIOD = 0.391878638976777*psr_utils.SECPERDAY # (in s)
ECCENTRICITY = 3.88136366443311e-05
PB_UNC = 8.2875e-11*psr_utils.SECPERDAY # Uncertainty on oribtal period (in s)
TSPAN = 667.203*psr_utils.SECPERDAY # Total span of timing solution (in s)

NSIG_DETECT = 3

def pbdot(pulsar_mass, companion_mass, pb, ecc):
    """Given the pulsar mass, companion mass, orbital period and eccenricity
        return the expected value of Pb-dot.

        Inputs:
            pulsar_mass: The mass of the pulsar (in Msun)
            companion_mass: The mass of the companion (in Msun)
            pb: The binary system's orbital period (in s)
            ecc: The binary system's eccentricity.

        Output:
            pbdot: The value of Pb-dot (in s/s).
    """
    f = lambda e: ((1 + (73.0/24)*e**2 + (37.0/96.0)*e**4)/((1-e**2)**(3.5)))
    pbdot = -((192*np.pi)/5.0)*((TSUN*2*np.pi)/pb)**(5.0/3.0)*f(ecc) * \
                ((pulsar_mass*companion_mass)/(pulsar_mass+companion_mass)**(1.0/3.0))
    return pbdot


def main():
    nump = 1000 # Number of pulsar masses to use
    numc = 1000 # Number of companion masses to use
    pulsar_masses = np.linspace(MP_MIN, MP_MAX, nump)
    comp_masses = np.linspace(MC_MIN, MC_MAX, numc)
   
    pbdots = np.empty((numc, nump))
    mp, mc = np.meshgrid(pulsar_masses, comp_masses)
    pbdots = pbdot(mp, mc, ORBITAL_PERIOD, ECCENTRICITY) # (in s/s)

    delta_pbdots = np.abs(pbdots*TSPAN) # (in s)
    numsigmas = delta_pbdots/PB_UNC
    numsigmas[numsigmas > NSIG_DETECT] = np.nan

    # Pb-dot is negative, but timespan to detect it is positive - take absolute value
    tspans_needed = np.abs(NSIG_DETECT*PB_UNC/pbdots) # (in s)

    # Set time spans where we should already see evidence for Pb-dot to NaN
    tspans_needed[tspans_needed < TSPAN] = np.nan

    # Plot
    fig = plt.figure(figsize=(8.5, 11))
    ax = plt.axes()
    plt.imshow(tspans_needed/psr_utils.SECPERDAY, origin='lower', aspect='auto', \
    #plt.imshow(numsigmas, origin='lower', aspect='auto', \
                extent=(pulsar_masses.min(), pulsar_masses.max(), \
                        comp_masses.min(), comp_masses.max()))
    cb = plt.colorbar(format=matplotlib.ticker.FuncFormatter(lambda val, ii: r"%d" % val))
    cb.set_label(r"Time span needed to detect $\.P_b$ (with $\sigma$=%d; days)" % NSIG_DETECT)
    #cb.set_label(r"Strength of Pb-dot detection today (Timing solution span: %d days; in $\sigma$)" % \
    #                (TSPAN/psr_utils.SECPERDAY))

    plt.axis([MP_MIN, MP_MAX, MC_MIN, MC_MAX]) # Set axes limits

    plt.xlabel(r"Pulsar Mass $M_p (M_\odot)$")
    plt.ylabel(r"Companion Mass $M_c (M_\odot)$")
    
    display = lambda mp, mc: r"Mp=%g, Mc=%g (tspan=%d days, Pb-dot=%.3g s/s)" % \
                (mp, mc, np.abs(NSIG_DETECT*PB_UNC/pbdot(mp,mc,ORBITAL_PERIOD, ECCENTRICITY)/psr_utils.SECPERDAY), \
                                    pbdot(mp,mc,ORBITAL_PERIOD, ECCENTRICITY))
    ax.format_coord = display
    fig.canvas.mpl_connect("key_press_event", \
                        lambda e: ((e.key in ('q', 'Q')) and plt.close(fig)))
    plt.show()


if __name__=="__main__":
    main()
