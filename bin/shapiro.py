#!/usr/bin/env python

"""
shapiro.py

A script to determine the size of shapiro delay for a low-eccentricity
orbit as a function of pulsar mass and companion mass.

Use equations 8.41, 8.50, 8.51 from Lorimer and Kramer's Pulsar Handbook.

Patrick Lazarus
Nov. 22, 2010
"""
import warnings

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker

import psr_utils

# Define a few constants
warnings.warn("Using hardcoded values for TRES, MASS_FUNC and PHI!")
TRES = 50e-6 # RMS residuals of fit in seconds
MASS_FUNC = 0.1531843160
PHI = np.pi/2.0


def sini(pulsar_mass, comp_mass, mass_func=MASS_FUNC):
    """Given a pulsar mass, companion mass and mass function return
        the corresponding sini.

        Inputs:
            pulsar_mass: Pulsar's mass (in solar mass units).
            comp_mass: Companion's mass (in solar mass units).
            mass_func: Mass function (in solar mass units).
    """
    sini = ((mass_func*(pulsar_mass+comp_mass)**2.0)**(1.0/3.0))/comp_mass
    return sini


def shapiro_delay(pulsar_mass, comp_mass, mass_func=MASS_FUNC, phi=PHI):
    """Given a pulsar mass and companion mass compute and return
        the shapiro delay.

        Inputs:
            pulsar_mass: Pulsar's mass (in solar mass units).
            comp_mass: Companion's mass (in solar mass units).
            mass_func: Mass function (in solar mass units).
            phi: Orbital phase measured from the ascending node.
    """
    warnings.warn("Assuming a low-eccentricity orbit!")
    range = psr_utils.Tsun*comp_mass
    shape = sini(pulsar_mass, comp_mass, mass_func)
    delay = -2*range*np.log(1-shape*np.sin(phi))
    return delay


def measurable_shapiro_delay(pulsar_mass, comp_mass, mass_func=MASS_FUNC, \
                                phi=PHI):
    """Given a pulsar mass and companion mass compute and return
        the measurable part of the shapiro delay (harmonics 3 and up).

        Inputs:
            pulsar_mass: Pulsar's mass (in solar mass units).
            comp_mass: Companion's mass (in solar mass units).
            mass_func: Mass function (in solar mass units).
            phi: Orbital phase measured from the ascending node.
    """
    warnings.warn("Assuming a low-eccentricity orbit!")
    range = psr_utils.Tsun*comp_mass
    shape = sini(pulsar_mass, comp_mass, mass_func)
    cbar = np.sqrt(1-shape**2)
    sigma = shape/(1+cbar) # Eq. 12 of Freire and Wex, 2010
    h3 = range*sigma**3 # Eq. 20 of F&W, 2010
    # The following is the measurable part of the Shapiro delay
    # Eq. 19 of F&W, 2010
    # measurable_delay = -4*h3*(np.sin(3*phi)/3.0 - sigma*np.cos(4*phi)/4.0 - \
    #                            sigma**2*np.sin(5*phi)/5.0 + \
    #                            sigma**3*np.cos(6*phi)/6.0)
    # use the exact equation: Eq. 28 of F&W, 2010
    measurable_delay = -2*h3*(np.log(1+sigma**2-2*sigma*np.sin(phi))/sigma**3 + \
                        2*np.sin(phi)/sigma**2 - np.cos(2*phi)/sigma)
    return measurable_delay


def main():
    nump = 1000 # Number of pulsar masses to use
    numc = 1000 # Number of companion masses to use
    pulsar_masses = np.linspace(1.2, 3.0, nump)
    comp_masses = np.linspace(0.9, 3.0, numc)
   
    delays = np.empty((numc, nump))
    inclination = np.empty((numc, nump))
    mp, mc = np.meshgrid(pulsar_masses, comp_masses)
    delays = measurable_shapiro_delay(mp, mc)
    inclination = np.arcsin(sini(mp, mc))*psr_utils.RADTODEG
    
    # Set delays that should already be detectable to an invalid number
    delays[delays > TRES] = np.nan 
    inclination[np.isnan(inclination)] = 91
    fig = plt.figure(figsize=(8.5, 11))
    # fig = plt.figure(figsize=(8,6.5))
    ax = plt.axes([0.1, 0.35, 0.85, 0.6])
    plt.imshow(np.log10(delays), origin='lower', aspect='auto', \
                extent=(pulsar_masses.min(), pulsar_masses.max(), \
                        comp_masses.min(), comp_masses.max()))
    cb = plt.colorbar(format=matplotlib.ticker.FuncFormatter(lambda val, ii: r"%4.1f" % (10**(6+val))))
    cb.set_label("Shapiro Delay Signal ($\mu s$)")
    contours = plt.contour(inclination, [90, 60, 45, 30], origin='lower', colors='k', \
                extent=(pulsar_masses.min(), pulsar_masses.max(), \
                        comp_masses.min(), comp_masses.max()))
    plt.clabel(contours, fmt=r"%d$^\circ$")
    
    hilite, = plt.plot(0,0, 'kx', mew='2')
    hilite.set_visible(False)
    
    plt.axis([1.2, 3.0, 0.9, 3.0]) # Set axes limits

    plt.xlabel(r"Pulsar Mass $M_p (M_\odot)$")
    plt.ylabel(r"Companion Mass $M_c (M_\odot)$")
    
    display = lambda mp, mc: "Mp=%g, Mc=%g (sini=%g, delay=%g us)" % \
                                (mp, mc, sini(mp,mc), \
                                    measurable_shapiro_delay(mp,mc))
    ax.format_coord = display

    ax2 = plt.axes([0.1, 0.05, 0.85, 0.25])
    ax2.format_coord = lambda phi, d: "Delay=%g us, Phi=%g" % \
                                (d, phi)
    line, = plt.plot(np.linspace(0,1,1000), np.zeros(1000), 'k-')
    plt.xlabel("Orbital Phase")
    plt.ylabel("Shapiro Delay ($\mu$s)")
    line.set_visible(False) 
    fig.canvas.mpl_connect("button_press_event", \
                            lambda e: show_shapiro_delay(e, ax, hilite, line))
    fig.canvas.mpl_connect("key_press_event", \
                        lambda e: ((e.key in ('q', 'Q')) and plt.close(fig)))
    plt.show()


def update_shapiro_signal(ax, hilite, line, pulsar_mass, comp_mass, mass_func=MASS_FUNC):
    """Update y-data in axes, 'ax', with the Shapiro delay signal
        for the given masses and mass function.

        Inputs:
            pulsar_mass: Pulsar's mass (in solar mass units).
            comp_mass: Companion's mass (in solar mass units).
            mass_func: Mass function (in solar mass units).
    """
    hilite.set_xdata(pulsar_mass)
    hilite.set_ydata(comp_mass)
    hilite.set_visible(True)

    phis = line.get_xdata()
    delays = measurable_shapiro_delay(pulsar_mass, comp_mass, \
                                            mass_func, phi=phis*np.pi*2)
    line.set_ydata(delays*1e6)
    line.axes.relim()
    line.axes.autoscale_view(scalex=False, scaley=True)
    line.set_visible(True)
    plt.draw()


def show_shapiro_delay(event, ax, hilite, line, mass_func=MASS_FUNC):
    if event.button==1 and event.inaxes==ax:
        update_shapiro_signal(ax, hilite, line, event.xdata, event.ydata, mass_func)

if __name__=="__main__":
    main()
