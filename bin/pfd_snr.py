#!/usr/bin/env python
"""
Manually determine on-pulse region for the profile in a pfd file
    and calculate the SNR.

    Patrick Lazarus, May 14, 2013 (taken out of gridding.py).
"""

import sys
import argparse
import warnings
import os.path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.integrate import trapz
from scipy.optimize import leastsq
import prepfold
import psr_utils
import injectpsr

from pypulsar.utils.astro import protractor
from pypulsar.utils.astro import sextant
from pypulsar.utils import estimate_snr
from pypulsar.utils import skytemp
debug = 1


class OnPulseError(Exception):
    pass


def transform(data, rot, scale, dc):
    nrot = int(np.round(rot*len(data)))
    #print "Rotating model by %d bins" % nrot
    rotated = np.asarray(psr_utils.rotate(data, nrot))
    return rotated*scale + dc


def get_rotation(profdata, modeldata, scale=1, dc=0):
    bestrot = 0
    bestrms = np.inf
    for rot in np.linspace(0, 1, len(profdata), endpoint=False):
        transformed = transform(modeldata, rot, scale, dc)
        resids = profdata - transformed
        rms = np.sqrt(np.mean(resids**2))
        if rms < bestrms:
            bestrms = rms
            bestrot = rot
    return bestrot


def get_resids(profdata, modeldata, scale=1, dc=0):
    if len(profdata) != len(modeldata):
        raise ValueError("Model and profile have different number " 
                         "of data points (%d and %d respectively)" % 
                         (len(modeldata), len(profdata)))
    bestrot = get_rotation(profdata, modeldata, scale)
    transformed = transform(modeldata, bestrot, scale, dc)
    resids = profdata - transformed
    return resids


def find_scale_and_phase(profdata, modeldata):
    to_optimize = lambda (scale, dc): get_resids(profdata, modeldata, scale, dc)
    init_params = [1, 0] # Initial multiplicative scale factor
    fit = leastsq(to_optimize, init_params)
    return fit


def read_gaussfitfile(gaussfitfile, proflen):
    """
    read_gaussfitfile(gaussfitfile, proflen):
        Read a Gaussian-fit file as created by the output of pygaussfit.py.
            The input parameters are the name of the file and the number of
            bins to include in the resulting template file.  A numpy array
            of that length is returned.

    Originally from PRESTO's psr_utils.py
    """
    phass = []
    ampls = []
    fwhms = []
    const = 0
    for line in open(gaussfitfile):
        if line.lstrip().startswith("phas"):
            phass.append(float(line.split()[2]))
        if line.lstrip().startswith("ampl"):
            ampls.append(float(line.split()[2]))
        if line.lstrip().startswith("fwhm"):
            fwhms.append(float(line.split()[2]))
        if line.lstrip().startswith("const"):
            const = float(line.split()[2])
    if not (len(phass) == len(ampls) == len(fwhms)):
        print "Number of phases, amplitudes, and FWHMs are not the same in '%s'!"%gaussfitfile
        return 0.0
    phass = np.asarray(phass)
    ampls = np.asarray(ampls)
    fwhms = np.asarray(fwhms)
    
    gauss_data = np.zeros((len(ampls), proflen))
    for ii in range(len(ampls)):
        data = ampls[ii]*psr_utils.gaussian_profile(proflen, phass[ii], fwhms[ii])
        dc = np.min(data)
        const += dc
        gauss_data[ii] = data - dc
    return gauss_data, const

class ObservationWithModel:
    """ ObservationWithModel object
    """
    def __init__(self, pfd, modelfn, sefd=None, verbose=True, bw=None):
        """Return an observation object for the given pfd file.
    
        Inputs: 
            pfd: a pfd object
            modelfn: the name of a file containing profile component parameters
            sefd: the system-equivalent flux density of the observation (in Jy)
            verbose: if True, print extra information (Default: True)
            bw: The bandwidth of the observation in MHz. (Default: determine from pfd file)

        Output:
            obs: The ObservationWithModel object            
        """
        self.sefd = sefd
        self.p = pfd 
        self.fn = self.p.pfd_filename
        self.snr = None
        self.smean = None
        self.verbose = verbose
        if bw is not None:
            self.bw = bw
        else:
            self.bw = self.p.chan_wid*self.p.numchan
        self.notes = []
        self.on_pulse = set()
        self.ignore_regions = []

        # Read model
        self.modelfn = modelfn
        self.modelparams = injectpsr.parse_model_file(self.modelfn)
        self.modelcomps = injectpsr.create_vonmises_components(self.modelparams)

        to_ignore = self.read_onpulse_file(self.modelfn+".on")
        self.p.dedisperse(doppler=True)
        self.p.adjust_period()
        if self.p.bestprof:
            prof = self.p.bestprof.profile
        else:
            prof = self.p.sumprof
        self.proflen = len(prof)
        self.nbin = len(prof)
        self.prof = np.asarray(prof)
        
        self.region_start = None
        self.region_start_line = None
        self.regions = []
       
        # Create components
        binphase = 1.0/self.nbin
        # Central phase of each bin
        phases = np.linspace(0, 1, self.nbin, endpoint=False) + 0.5*binphase
        modeldata = np.zeros(len(self.prof))
        for vm in self.modelcomps:
            compdata = vm(phases)
            modeldata += compdata

        # Fit individual analytic components to profile
        params, fitcode = find_scale_and_phase(self.prof, modeldata)
        scale, dc = params
        rot = get_rotation(self.prof, modeldata, scale)
        
        # Plot
        self.fig = plt.gcf()
        self.profax = plt.subplot(3,1,1)
        plt.plot(self.prof-dc, 'k-', drawstyle='steps-post')
        plt.title("Profile")
        
        # Individual components
        self.modelax = plt.subplot(3,1,2, sharey=self.profax)
        self.compartists = []
        self.comp_sums = []
        self.comp_maxes = []
        for ii, vm in enumerate(self.modelcomps):
            compdata = vm(phases)
            transformed = transform(compdata, rot, scale, 0)
            lw = 1+int(ii in self.on_pulse)
            self.compartists.append(plt.plot(transformed, ls='-', 
                                             drawstyle='steps', 
                                             lw=lw, picker=True, 
                                             label="Comp. #%d" % ii)[0])
            self.comp_sums.append(np.sum(transformed))
            self.comp_maxes.append(np.max(transformed))
        plt.legend(loc='best', prop=dict(size='x-small'))
        modeldata = transform(modeldata, rot, scale, dc)
        
        # Plot residuals
        self.residax = plt.subplot(3,1,3)
        self.residuals = self.prof - modeldata
        plt.plot(self.residuals, c='#444444', ls='-', drawstyle='steps-post')
        plt.axhline(0, c='k', ls='--')
        for xlo, xhi in to_ignore:
            self.regions.append(self.residax.axvspan(xlo, xhi, picker=True,
                                            facecolor='r', alpha=0.5))
            

        # Set up triggers
        self.cid_pick = self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.cid_keypress = self.fig.canvas.mpl_connect('key_press_event', \
                                                            self.keypress)
        self.cid_mousepress = self.fig.canvas.mpl_connect('button_press_event', \
                                                            self.mousepress)
    
    def read_onpulse_file(self, fn):
        # Read on-pulse components
        to_ignore = []
        if os.path.exists(fn):
            self.read_onpulse_from_file = True
            with open(fn, 'r') as ff:
                for line in ff:
                    line, sep, comment = line.partition('#')
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith("ignore:"):
                        xlo, xhi = line.split()[1:3]
                        to_ignore.append((int(xlo), int(xhi)))
                    else:
                        self.on_pulse.add(int(line))
        else:
            self.read_onpulse_from_file = False
        self.read_onpulse_from_file = True
        return to_ignore

    def onpick(self, event):
        if (event.mouseevent.inaxes == self.modelax) and \
                    (event.mouseevent.button == 1):
            ind = self.compartists.index(event.artist)
            if ind in self.on_pulse:
                print "Component %d un-selected" % ind
                self.on_pulse.remove(ind)
                event.artist.set_linewidth(1)
            else:
                print "Component %d selected as on-pulse" % ind
                self.on_pulse.add(ind)
                event.artist.set_linewidth(2)
            self.fig.canvas.draw()
        elif (event.mouseevent.inaxes == self.residax) and \
                    (event.mouseevent.button == 3): 
            ind = self.regions.index(event.artist)
            self.regions.pop(ind)
            event.artist.remove()
            self.fig.canvas.draw()

    def calc_snr(self):
        if len(self.on_pulse) < 1:
            warnings.warn("No on-pulse region selected!")
            return
        area = 0
        profmax = 0
        if self.on_pulse or self.regions:
            with open(self.modelfn+".on", 'w') as ff:
                ff.write("#On-pulse components:\n")
                for ionpulse in self.on_pulse:
                    ff.write("%d\n" % ionpulse)
                for polygon in self.regions:
                    xlo, xhi = polygon.get_xy()[0:3:2,0]
                    ff.write("ignore: %d %d\n" % (xlo, xhi))
        for ionpulse in self.on_pulse:
            area += self.comp_sums[ionpulse]
            profmax = max(profmax, self.comp_maxes[ionpulse])
        
        # Correct standard deviation for correlations between bins
        #data_avg, data_var = self.p.stats.sum(axis=1).mean(axis=0)[1:3]
        #nbin_eff = self.proflen*self.p.DOF_corr()
        #std = np.sqrt(data_var*self.p.Nfolded/nbin_eff)
       
        iignore = np.zeros_like(self.residuals, dtype=bool)
        for polygon in self.regions:
            xlo, xhi = polygon.get_xy()[0:3:2,0]
            iignore[xlo:xhi] = 1
        std = np.std(self.residuals[~iignore])
        std /= self.p.DOF_corr()
       
        # Calculate S/N using eq. 7.1 from Lorimer and Kramer
        self.weq = area/profmax
        self.snr = area/std/np.sqrt(self.weq)
        if self.verbose:
            if debug:
                print "Equivalent width (bins):", self.weq
                print "Std-dev correction factor:", self.p.DOF_corr()
                print "Std-dev corrected for correlations between phase bins:", std
                print "Integral under the selected pulse components:", \
                        area
            print "SNR:", self.snr

        if self.sefd is not None:
            npol = 2  # prepfold files only contain total-intensity
                      # (i.e. both polarisations summed)
            self.smean = self.snr*self.sefd/np.sqrt(npol*self.p.T*self.bw)*np.sqrt(self.weq/(len(self.prof)-self.weq))
            if self.verbose:
                print "Mean flux density (mJy):", self.smean
    
    def mousepress(self, event):
        if (event.inaxes == self.residax) and event.button==1:
            self.eventpress = event
            if self.region_start is None:
                # Starting a new ignore region
                xx =  int(event.xdata+0.5)
                self.region_start = xx
                self.region_start_line = plt.axvline(xx, c='k', ls='--')
            else:
                # Selected an on-pulse region
                if self.region_start > event.xdata:
                    xlo = 0
                    xhi = int(event.xdata+0.5)
                    self.regions.append(self.residax.axvspan(xlo, xhi, picker=True,
                                                    facecolor='r', alpha=0.5))
                    xlo = int(self.region_start+0.5)
                    xhi = self.nbin
                    self.regions.append(self.residax.axvspan(xlo, xhi, picker=True,
                                                    facecolor='r', alpha=0.5))
                else:
                    xlo = int(self.region_start+0.5)
                    xhi = int(event.xdata+0.5)
                    self.regions.append(self.residax.axvspan(xlo, xhi, picker=True,
                                                    facecolor='r', alpha=0.5))
                
                # Remove line
                self.region_start_line.remove()
                self.region_start = None
                self.region_start_line = None

            self.fig.canvas.draw()
    
    def keypress(self, event):
        if event.key == ' ':
            self.calc_snr()
        elif event.key in ('q', 'Q'):
            self.calc_snr()
            plt.close(self.fig)
        elif event.key in ('r', 'R'):
            self.notes.append("RFI!")
        elif event.key in ('n', 'N'):
            self.notes.append("No detection!")


class ObservationWithGauss:
    """ ObservationWithGauss object
    """
    def __init__(self, pfd, gaussfn, sefd=None, verbose=True, bw=None):
        """Return an observation object for the given pfd file.
    
        Inputs: 
            pfd: a pfd object
            gaussfn: the name of a file containing profile component parameters
            sefd: the system-equivalent flux density of the observation (in Jy)
            verbose: if True, print extra information (Default: True)
            bw: The bandwidth of the observation in MHz. (Default: determine from pfd file)

        Output:
            obs: The ObservationWithModel object            
        """
        self.sefd = sefd
        self.p = pfd 
        self.fn = self.p.pfd_filename
        self.snr = None
        self.smean = None
        self.verbose = verbose
        if bw is not None:
            self.bw = bw
        else:
            self.bw = self.p.chan_wid*self.p.numchan
        self.notes = []
        self.on_pulse = set()
        self.ignore_regions = []

        self.p.dedisperse(doppler=True)
        self.p.adjust_period()
        if self.p.bestprof:
            prof = self.p.bestprof.profile
        else:
            prof = self.p.sumprof
        self.proflen = len(prof)
        self.nbin = len(prof)
        self.prof = np.asarray(prof)
        
        self.region_start = None
        self.region_start_line = None
        self.regions = []
       
        # Read gaussians file
        self.gaussfn = gaussfn
        self.gauss_data, dc = read_gaussfitfile(self.gaussfn, self.nbin)
        modeldata = np.sum(self.gauss_data, axis=0)+dc
        to_ignore = self.read_onpulse_file(self.gaussfn+".on")

        # Plot
        self.fig = plt.gcf()
        self.profax = plt.subplot(3,1,1)
        plt.plot(self.prof-dc, 'k-', drawstyle='steps-post')
        plt.title("Profile")
        
        # Individual components
        self.modelax = plt.subplot(3,1,2, sharey=self.profax)
        self.compartists = []
        self.comp_sums = []
        self.comp_maxes = []
        for ii, gauss in enumerate(self.gauss_data):
            lw = 1+int(ii in self.on_pulse)
            self.compartists.append(plt.plot(gauss, ls='-', 
                                             drawstyle='steps-pre', 
                                             lw=lw, picker=True, 
                                             label="Comp. #%d" % ii)[0])
            self.comp_sums.append(np.sum(gauss))
            self.comp_maxes.append(np.max(gauss))
        plt.legend(loc='best', prop=dict(size='x-small'))
        
        # Plot residuals
        self.residax = plt.subplot(3,1,3)
        self.residuals = self.prof - modeldata
        plt.plot(self.residuals, c='#444444', ls='-', drawstyle='steps-post')
        plt.axhline(0, c='k', ls='--')
        for xlo, xhi in to_ignore:
            self.regions.append(self.residax.axvspan(xlo, xhi, picker=True,
                                            facecolor='r', alpha=0.5))
            
        # Set up triggers
        self.cid_pick = self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.cid_keypress = self.fig.canvas.mpl_connect('key_press_event', \
                                                            self.keypress)
        self.cid_mousepress = self.fig.canvas.mpl_connect('button_press_event', \
                                                            self.mousepress)
    
    def read_onpulse_file(self, fn):
        # Read on-pulse components
        to_ignore = []
        if os.path.exists(fn):
            self.read_onpulse_from_file = True
            with open(fn, 'r') as ff:
                for line in ff:
                    line, sep, comment = line.partition('#')
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith("ignore:"):
                        xlo, xhi = line.split()[1:3]
                        to_ignore.append((int(xlo), int(xhi)))
                    else:
                        self.on_pulse.add(int(line))
        else:
            self.read_onpulse_from_file = False
        self.read_onpulse_from_file = True
        return to_ignore

    def onpick(self, event):
        if (event.mouseevent.inaxes == self.modelax) and \
                    (event.mouseevent.button == 1):
            ind = self.compartists.index(event.artist)
            if ind in self.on_pulse:
                print "Component %d un-selected" % ind
                self.on_pulse.remove(ind)
                event.artist.set_linewidth(1)
            else:
                print "Component %d selected as on-pulse" % ind
                self.on_pulse.add(ind)
                event.artist.set_linewidth(2)
            self.fig.canvas.draw()
        elif (event.mouseevent.inaxes == self.residax) and \
                    (event.mouseevent.button == 3): 
            ind = self.regions.index(event.artist)
            self.regions.pop(ind)
            event.artist.remove()
            self.fig.canvas.draw()

    def calc_snr(self):
        if len(self.on_pulse) < 1:
            warnings.warn("No on-pulse region selected!")
            return
        area = 0
        profmax = 0
        if self.on_pulse or self.regions:
            with open(self.gaussfn+".on", 'w') as ff:
                ff.write("#On-pulse components:\n")
                for ionpulse in self.on_pulse:
                    ff.write("%d\n" % ionpulse)
                for polygon in self.regions:
                    xlo, xhi = polygon.get_xy()[0:3:2,0]
                    ff.write("ignore: %d %d\n" % (xlo, xhi))
        for ionpulse in self.on_pulse:
            area += self.comp_sums[ionpulse]
            profmax = max(profmax, self.comp_maxes[ionpulse])
        
        # Correct standard deviation for correlations between bins
        #data_avg, data_var = self.p.stats.sum(axis=1).mean(axis=0)[1:3]
        #nbin_eff = self.proflen*self.p.DOF_corr()
        #std = np.sqrt(data_var*self.p.Nfolded/nbin_eff)
       
        iignore = np.zeros_like(self.residuals, dtype=bool)
        for polygon in self.regions:
            xlo, xhi = polygon.get_xy()[0:3:2,0]
            iignore[xlo:xhi] = 1
        std = np.std(self.residuals[~iignore])
        std /= self.p.DOF_corr()
       
        # Calculate S/N using eq. 7.1 from Lorimer and Kramer
        self.weq = area/profmax
        self.snr = area/std/np.sqrt(self.weq)
        if self.verbose:
            if debug:
                print "Equivalent width (bins):", self.weq
                print "Std-dev correction factor:", self.p.DOF_corr()
                print "Std-dev corrected for correlations between phase bins:", std
                print "Integral under the selected pulse components:", \
                        area
            print "SNR:", self.snr

        if self.sefd is not None:
            npol = 2  # prepfold files only contain total-intensity
                      # (i.e. both polarisations summed)
            self.smean = self.snr*self.sefd/np.sqrt(npol*self.p.T*self.bw)*np.sqrt(self.weq/(len(self.prof)-self.weq))
            if self.verbose:
                print "Mean flux density (mJy):", self.smean
    
    def mousepress(self, event):
        if (event.inaxes == self.residax) and event.button==1:
            self.eventpress = event
            if self.region_start is None:
                # Starting a new ignore region
                xx =  int(event.xdata+0.5)
                self.region_start = xx
                self.region_start_line = plt.axvline(xx, c='k', ls='--')
            else:
                # Selected an on-pulse region
                if self.region_start > event.xdata:
                    xlo = 0
                    xhi = int(event.xdata+0.5)
                    self.regions.append(self.residax.axvspan(xlo, xhi, picker=True,
                                                    facecolor='r', alpha=0.5))
                    xlo = int(self.region_start+0.5)
                    xhi = self.nbin
                    self.regions.append(self.residax.axvspan(xlo, xhi, picker=True,
                                                    facecolor='r', alpha=0.5))
                else:
                    xlo = int(self.region_start+0.5)
                    xhi = int(event.xdata+0.5)
                    self.regions.append(self.residax.axvspan(xlo, xhi, picker=True,
                                                    facecolor='r', alpha=0.5))
                
                # Remove line
                self.region_start_line.remove()
                self.region_start = None
                self.region_start_line = None

            self.fig.canvas.draw()
    
    def keypress(self, event):
        if event.key == ' ':
            self.calc_snr()
        elif event.key in ('q', 'Q'):
            self.calc_snr()
            plt.close(self.fig)
        elif event.key in ('r', 'R'):
            self.notes.append("RFI!")
        elif event.key in ('n', 'N'):
            self.notes.append("No detection!")


class Observation:
    """ Observation object
    """
    def __init__(self, pfd, sefd=None, verbose=True):
        """Return an observation object for the given pfd file.
    
        Inputs: 
            pfd: a pfd object
            sefd: the system-equivalent flux density of the observation (in Jy)
            verbose: if True, print extra information (Default: True)

        Output:
            obs: The Observation object            
        """
        self.sefd = sefd
        self.p = pfd 
        self.fn = self.p.pfd_filename
        self.snr = None
        self.smean = None
        self.verbose = verbose
        self.notes = []
        
        self.p.dedisperse(doppler=True)
        self.p.adjust_period()
        if self.p.bestprof:
            prof = self.p.bestprof.profile
        else:
            prof = self.p.sumprof
        self.proflen = len(prof)
        self.nbin = len(prof)
        imax = np.argmax(prof)
        self.nrot = (imax-len(prof)/2) % len(prof)
        if self.verbose:
            print "Profile maximum at bin %d. Rotating by %d bins." % (imax, self.nrot)
        self.prof = np.asarray(psr_utils.rotate(prof, self.nrot))
        
        self.region_start = None
        self.region_start_line = None
        self.regions = []
        
        # Plot
        self.fig = plt.gcf()
        self.ax = plt.gca()
        plt.plot(self.prof, 'k-', drawstyle='steps-post')
        
        # Set up triggers
        self.cid_mouseclick = self.fig.canvas.mpl_connect('button_press_event', \
                                                            self.mousepress)
        self.cid_pick = self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.cid_keypress = self.fig.canvas.mpl_connect('key_press_event', \
                                                            self.keypress)

    def mousepress(self, event):
        if event.inaxes and event.button==1:
            self.eventpress = event
            if self.region_start is None:
                # Starting a new on-pulse region
                xx =  int(event.xdata+0.5)
                self.region_start = xx
                self.region_start_line = plt.axvline(xx, c='k', ls='--')
            else:
                # Selected an on-pulse region
                if self.region_start > event.xdata:
                    xlo = 0
                    xhi = int(event.xdata+0.5)
                    self.regions.append(plt.axvspan(xlo, xhi, picker=True,
                                                    facecolor='b', alpha=0.5))
                    xlo = int(self.region_start+0.5)
                    xhi = self.nbin
                    self.regions.append(plt.axvspan(xlo, xhi, picker=True,
                                                    facecolor='b', alpha=0.5))
                else:
                    xlo = int(self.region_start+0.5)
                    xhi = int(event.xdata+0.5)
                    self.regions.append(plt.axvspan(xlo, xhi, picker=True,
                                                    facecolor='b', alpha=0.5))
                
                # Remove line
                self.region_start_line.remove()
                self.region_start = None
                self.region_start_line = None

            self.fig.canvas.draw()

    def onpick(self, event):
        if event.mouseevent.button==3:
            ind = self.regions.index(event.artist)
            self.regions.pop(ind)
            event.artist.remove()
            self.fig.canvas.draw()

    def calc_snr(self):
        ionpulse = np.zeros_like(self.prof, dtype=bool)
        for polygon in self.regions:
            xlo, xhi = polygon.get_xy()[0:3:2,0]
            ionpulse[xlo:xhi] = 1

        nbins_selected = ionpulse.sum()
        if nbins_selected < 1:
            warnings.warn("No on-pulse region selected!")
            return
        # Correct standard deviation for correlations between bins
        data_avg, data_var = self.p.stats.sum(axis=1).mean(axis=0)[1:3]
        #print data_var
        nbin_eff = self.proflen*self.p.DOF_corr()
        std = np.sqrt(data_var*self.p.Nfolded/nbin_eff)
       
        # Calculate S/N using eq. 7.1 from Lorimer and Kramer
        offpulse = self.prof[~ionpulse]

        mean = offpulse.mean()
        scaled = self.prof-mean
        area = np.sum(scaled[ionpulse]) 
        profmax = np.max(scaled[ionpulse])
        self.weq = area/profmax
        self.snr = area/std/np.sqrt(self.weq)
        if self.verbose:
            print "Number of bins selected: %d (%f phase)" % \
                    (nbins_selected, nbins_selected/float(len(self.prof)))
            if debug:
                print "Equivalent width (bins):", self.weq
                print "Std-dev corrected for correlations between phase bins:", std
                print "Off-pulse mean:", mean
                print "Integral under the mean-subtracted on-pulse region:", \
                        area
            print "SNR:", self.snr

        if self.sefd is not None:
            npol = 2  # prepfold files only contain total-intensity
                      # (i.e. both polarisations summed)
            bw = self.p.chan_wid*self.p.numchan
            self.smean = self.snr*self.sefd/np.sqrt(npol*self.p.T*bw)*np.sqrt(self.weq/(len(self.prof)-self.weq))
            if self.verbose:
                print "Mean flux density (mJy):", self.smean
    
    def keypress(self, event):
        if event.key == ' ':
            self.calc_snr()
        elif event.key in ('q', 'Q'):
            self.calc_snr()
            plt.close(self.fig)
        elif event.key in ('r', 'R'):
            self.notes.append("RFI!")
        elif event.key in ('n', 'N'):
            self.notes.append("No detection!")


def main():
    for pfdfn in args.files:
        print pfdfn
        pfd = prepfold.pfd(pfdfn)
        if args.sefd is not None:
            sefd = args.sefd
        elif args.gain is not None and args.tsys is not None:
            fctr = 0.5*(pfd.hifreq+pfd.lofreq)
            glon, glat = sextant.equatorial_to_galactic(pfd.rastr, pfd.decstr, 
                                                        input='sexigesimal', 
                                                        output='deg')
            print "Galactic Coords: l=%g deg, b=%g deg" % (glon, glat)
            tsky = skytemp.get_skytemp(glon, glat, freq=fctr)[0]
            print "Sky temp at %g MHz: %g K" % (fctr, tsky)
            sefd = (args.tsys+tsky)/args.gain
        print sefd, args.fwhm, args.sep
        if (sefd is not None) and (args.fwhm is not None) and (args.sep is not None):
            factor = estimate_snr.airy_pattern(args.fwhm, args.sep)
            print "Pulsar is off-centre"
            print "Reducing SEFD by factor of %g (SEFD: %g->%g)" % (factor, sefd, sefd/factor)
            sefd /= factor
        if args.model_file is not None:
            obs = ObservationWithModel(pfd, args.model_file, sefd=sefd, verbose=True)
        elif args.gauss_file is not None:
            obs = ObservationWithGauss(pfd, args.gauss_file, sefd=sefd, verbose=True)
        else:
            obs = Observation(pfd, sefd=sefd, verbose=True)
        plt.show()
        print " ".join(obs.notes)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate SNR from .pfd files.")
    parser.add_argument('files', nargs='*', \
                        help="Files to find SNR for.")
    parser.add_argument('--on-pulse', dest='on_pulse', nargs=2, type=float, \
                        help="On-pulse region. Two values should be provided, " \
                            "the starting phase and the ending phase.")
    parser.add_argument('--sefd', dest='sefd', type=float, \
                        help="The SEFD (in Jy) of the observing system. "
                             "(i.e. Tsys/Gain) Sky temperature will not be factored in.")
    parser.add_argument('--tsys', dest='tsys', type=float, \
                        help="The temperature of the observing system in K "
                             "(not including sky temperature).")
    parser.add_argument('--gain', dest='gain', type=float, \
                        help="The gain of the observing system in K/Jy.")
    parser.add_argument('--sep', dest='sep', type=float, \
                        help="The angular separation, in arcmin, between the "
                             "beam centre and the pulsar. This reduces the "
                             "effective gain of the observation by assuming "
                             "an Airy disk. (The --fwhm option must also be provided.)")
    parser.add_argument('--fwhm', dest='fwhm', type=float, \
                        help="The FWHM of the beam's Airy disk pattern, in arcmin.")
    parser.add_argument('-m', '--model-file', dest='model_file', type=str, 
                        default=None,
                        help="A paas-created .m file containing parameters " 
                             "describing components fit to the profile.")
    parser.add_argument('-g', '--gaussian-file', dest='gauss_file', type=str, 
                        default=None,
                        help="A pygaussfit.py created gaussians file "
                             "containing parameters " 
                             "describing components fit to the profile.")
    args = parser.parse_args()

    if args.sefd is not None and (args.tsys is not None or args.gain is not None):
        raise ValueError("Gain and/or system temperature should not be " 
                         "provided if SEFD is given.")
    elif (args.tsys is not None and args.gain is None) or \
         (args.tsys is None and args.gain is not None):
        raise ValueError("Both gain and system temperature must be provided " 
                         "together.")
    main()
