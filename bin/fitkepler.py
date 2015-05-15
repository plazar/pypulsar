#!/usr/bin/env python

# A program to fit periods with a keplerian orbit.
#   (See http://www.scipy.org/Cookbook/FittingData "3.2 Fitting a 2D Gaussian")
#
#       Patrick Lazarus, March 20th, 2009

import sys
import optparse
import glob
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import prepfold
import psr_utils
import bestprof

count = 0

PARAMNAMES = ("Asini (lt-s)","Orbital period (days)","Pulsar spin period (s)",\
                "T0 (MJD)","Eccentricity","Longitude of periastron (rad)")

def read_pfds(fns):
    """
    Read pfd files and return arrays of period, period error and mjds
    
    Input: fns - a list of files that contain a list of pfd files to read
    """
    ps = []
    perrs = []
    mjds = []
    
    for fn in fns:
        file = open(fn, 'r')
        pfdfns = [x.strip() for x in file.readlines() if x.strip()[0] != '#']
        file.close()
        for pfdfn in pfdfns:
            p = prepfold.pfd(pfdfn)
            b = bestprof.bestprof("%s.bestprof" % pfdfn)
            if options.extendnum and (b.p1err_bary < b.p1_bary):
                # If extension is reqested and p-dot is at least a 4-sigma detection
                maxt = 300#np.sqrt((0.25*b.p0_bary)**2 - b.p0err_bary**2)/b.p1err_bary
                tstep = maxt/options.extendnum
                if options.debuglevel:
                    print "Bary Period (s): %g +/- %g" % (b.p0_bary, b.p0err_bary)
                    print "Bary P-dot (s/s): %g +/- %g" % (b.p1_bary, b.p1err_bary)
                    print "Bary Epoch (MJD): %.12f" % p.bepoch
                    print "    Period points (included phantom period measurements)"
                    print "    (T_max = %g (s), T_step = %g (s))" % (maxt, tstep)
                for ii in np.arange(-options.extendnum, options.extendnum+1):
                    ps.append(b.p0_bary+b.p1_bary*ii*tstep)
                    perrs.append(np.sqrt((b.p0err_bary)**2+(ii*tstep*b.p1err_bary)**2)*options.efac)
                    mjds.append(p.bepoch+ii*tstep/psr_utils.SECPERDAY)
                    print "  %.15f  %.10f   %.10f" % (mjds[-1], ps[-1]*1000, perrs[-1]*1000)
                    if options.debuglevel:
                        print "   P (s): %g +/- %g @ MJD=%.12f" % (ps[-1], perrs[-1], mjds[-1])
            else:
                ps.append(b.p0_bary)
                perrs.append(b.p0err_bary*options.efac)
                mjds.append(p.bepoch)
                print "  %.15f  %.10f   %.10f" % (mjds[-1], ps[-1]*1000, perrs[-1]*1000)
    return np.array(ps), np.array(perrs), np.array(mjds)

def read_textfile(fns):
    """
    Read text file which contains 3 columns: mjd, period, period error
    
    Input: filename
    """
    ps = []
    perrs = []
    mjds = []
    
    for fn in fns:
        file = open(fn, 'r')
        for ii, line in enumerate(file.readlines()):
            mjd, p, perr = line.strip().split()[:3]
            mjds.append(float(mjd))
            ps.append(float(p))
            perrs.append(float(perr)*options.efac)
    
    return np.array(ps), np.array(perrs), np.array(mjds)

def print_params(params):
    for name, val in zip(PARAMNAMES, params):
        print "\t%s: %r" % (name, val)

def kepler_function(asini, p_orb, p_psr, T0, ecc=0, peri=0):
    """
    Given amplitude, period and offset return a sine function.
    Asini is in light-seconds  (?)
    P_orb is in days.
    P_psr is in s.
    T0 is in mjd.
    Ecc is the eccentricity.
    Peri is the angle of periastron in radians.
    The function accepts one argument, mjd.
    """
    def func(mjd):
        global count
        count += 1
        if options.debuglevel > 1:
            print "Evaluation #:", count
            print "Params:"
            print_params((asini, p_orb, p_psr, T0, ecc, peri))
        if options.debuglevel > 2:
            print "\tMJD:", mjd
        p_orb_sec = p_orb*psr_utils.SECPERDAY # orbital period in seconds
        orb_freq_hz = psr_utils.TWOPI/p_orb_sec # Orbital angular frequency in hertz
        orb_freq = psr_utils.TWOPI/p_orb # Orbital angular frequency
        
        ma = between_zero_twopi(orb_freq*(mjd-T0)) # Mean anomaly
        E = between_zero_twopi(eccentric_anomaly(ecc, ma)) # Eccentric anomaly
        A = between_zero_twopi(2*np.arctan(np.sqrt((1+ecc)/(1-ecc))*np.tan(E/2.0))) # True anomaly
        
        velocity = orb_freq_hz*asini/np.sqrt(1-ecc**2)*(np.cos(peri+A)+ecc*np.cos(peri)) # in units of C
        observed_p_psr = p_psr*(1+velocity) # velocity in units of C
        
        if options.debuglevel > 2:
            print "Intermediate calculations:"
            print "\tp_orb_sec:", p_orb_sec
            print "\torb_freq_hz:", orb_freq_hz
            print "\torb_freq:", orb_freq
            print "\tma:", ma
            print "\tE:", E
            print "\tA:", A
            print "\tvelocity:", velocity
            print "\tobserved_p_psr:", observed_p_psr
        return observed_p_psr
    return func

def eccentric_anomaly(eccentricity, mean_anomaly):
    """
    eccentric_anomaly(eccentricity, mean_anomaly):
        Return the eccentric anomaly in radians, given a set of mean_anomalies
        in radians.
    
    Taken from Astronomical Algorithms by Jean Meeus
    """
    ma = between_zero_twopi(mean_anomaly)
    F = np.ones(mean_anomaly.size)
    F[ma > np.pi] = -1
    ma = np.where(ma > np.pi, psr_utils.TWOPI-ma, ma)
    D = np.pi/4.0
    ecc_anom = psr_utils.PIBYTWO
    for j in range(53):
        ma1 = ecc_anom-eccentricity*np.sin(ecc_anom)
        ecc_anom = ecc_anom + D*np.sign(ma-ma1)
        D = D/2.0
    return ecc_anom*F

def between_zero_twopi(rad):
    """
    return angles between 0 and 2pi
    """
    r = np.fmod(rad, psr_utils.TWOPI)
    return np.where(r < 0.0, r+psr_utils.TWOPI, r)

def min_comp_mass(Pb, x):
    """
    given Pb (days) and asini (lt-s) use scipy's newton-method root-solver
    to calculate minimum companion mass (assuming edge-on orbit and pulsar
    is 1.4 solar-masses.)
    """
    mp = 1.4 # pulsar mass (assume 1.4 Msol)
    i = psr_utils.PIBYTWO # inclination angle in rad (assume edge on for min companion mass)
    
    f = lambda mc: psr_utils.mass_funct(Pb*psr_utils.SECPERDAY,np.fabs(x))-psr_utils.mass_funct2(mp,mc,i)

    return opt.newton(f, 0.1) # 0.1 is initial guess

def fit(params, data):
    """
    returns the parameters found by a fit to data provided
    """
    ps, perrs, mjds = data
    errorfunction = lambda p, vals, errs: np.ravel((kepler_function(*p)(mjds) - vals)/errs)
    if options.debuglevel > 0:
         p, cov_x, infodict, mesg, success = opt.leastsq(errorfunction, params, maxfev=options.maxfev, full_output=1, args=(ps, perrs))
         print "Covariance matrix:\n", cov_x
    else:
         p, success = opt.leastsq(errorfunction, tuple(params), maxfev=options.maxfev, args=(ps, perrs))
    if success in [1,2,3,4]:
        # Fit is successful
        pass
    else:
        raise "Fit FAILED! Run program with 'debug' on to get more information"
    return p

def main():
    # Get filenames
    fns = glob.glob(options.filename)
    if options.usepfds==True:
        data = read_pfds(fns)
    else:
        print "reading from", fns
        data = read_textfile(fns)
    ps, perrs, mjds = data
    if options.debuglevel > 0:
        print "Data:"
        for z in zip(*data):
            print "\tPeriod:", z[0], "Period error:", z[1], "MJD:", z[2]
    if options.debuglevel > 0:
        print "initial parameters:"
        print_params(init_params) 
    params = init_params.copy() # copy of init_params
    print "Fitting %d data points" % len(mjds)
    result = fit(params, data)
    # Display results of fit
    print "Fit results:"
    print_params(result)
    print "\tMin companion mass: ", min_comp_mass(result[1], result[0])

    if len(predict_mjds) > 0:
        periods = kepler_function(*init_params)(predict_mjds)
        print "\nSpin period predictions:"
        for mjd, p in zip(predict_mjds,periods):
            print "\t%r: %r s" % (mjd, p)
    
    t_actual = np.linspace(mjds.min()-0.5*result[1], mjds.max()+0.5*result[1], mjds.ptp()*1000)
    t = t_actual - int(mjds.min())
    plt.figure(figsize=(11,8.5))
    ax = plt.subplot(2,1,1)
    if options.debuglevel > 0:
        # plot initial guess
        plt.plot(t, kepler_function(*init_params)(t_actual)-result[2], 'r--')
    plt.plot(t, kepler_function(*result)(t_actual)-result[2], 'k--')
    plt.axhline(0, ls=':', color='k')
    plt.errorbar(mjds-int(mjds.min()), ps-result[2], yerr=perrs, fmt='k.')
    plt.ylabel("Bary Period (s) - %f" % result[2])
    plt.xlabel("Epoch (MJD) - %d" % mjds.min())
    ax.ticklabel_format(style='plain', axes='both')
    plt.subplot(2,1,2, sharex=ax)
    kepler = kepler_function(*result)
    plt.errorbar(mjds-int(mjds.min()), ps-kepler(mjds), yerr=perrs, fmt='k.')
    plt.savefig(options.savefn, orientation='landscape', papertype='letter')
    plt.figure()
    dataphases = np.fmod(mjds-result[3], result[1])/result[1]
    dataphases[dataphases<0] = dataphases+1
    phasegrid = np.linspace(0,1,1000, endpoint=True)
    plt.plot(phasegrid, kepler(result[3]+phasegrid*result[1])-result[2], 'k--')
    if options.debuglevel > 0:
        # plot initial guess
        initphases = np.fmod(mjds-init_params[3], init_params[1])/init_params[1]
        initphases[initphases<0] = initphases+1
        plt.plot(phasegrid, kepler_function(*init_params)(init_params[3]+phasegrid*init_params[1])-init_params[2], 'r--')
        plt.errorbar(initphases, ps-init_params[2], yerr=perrs, fmt='r.')
    plt.errorbar(dataphases, ps-result[2], yerr=perrs, fmt='k.')
    plt.xlabel("Orbital Phase")
    plt.ylabel("Bary Period (s) - %f" % result[2])
    if options.plot:
        plt.show()

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('--use-pfds', dest='usepfds', action='store_true', help="Files provided contain a list of pfd files to read. (This is the default)", default=True)
    parser.add_option('--use-text', dest='usepfds', action='store_false', help="Files provided are text files containing 3 columns: mjd, period, period error. (Default is files contain list of pfd files to read)", default=True)
    parser.add_option('-f', '--file', dest='filename', type='string', help="File containing a list of pfd files to use, or a list of mjds and periods.", default='')
    parser.add_option('-d', '--debug', dest='debuglevel', type='int', help="Set the debug level. Higher level means more information will be printed to screen. (Default: No debug info", default=0)
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', help="Be verbose. Same as --debug=1. (Default: Not verbose)", default=False)
    parser.add_option('-x', '--extend', dest='extendnum', type='int', help="If barycentric p-dot is significant, estimate period just before and after measured value and include the estimated values in the fit. The argument should be the maximum number of estimated values to use on either side of the measured period. (Default: 0)", default=0)
    parser.add_option('-e', '--efac', dest='efac', type='float', help="Factor to multiply all uncertainties by. (Default: 1)", default=1)
    parser.add_option('--maxfev', dest='maxfev', type='int', help="Maximum number of function evaluations in fitting. (Default: 1000)", default=1000)
    parser.add_option('-s', '--savefn', dest='savefn', type='string', help="Filename to save plot to. (Default: 'fitkepler.tmp.ps')", default='fitkepler.tmp.ps')
    parser.add_option('-n', '--noplot', dest='plot', action='store_false', help="Do not show plot. Only make ps file. (Default: Show plot)", default=True)
    parser.add_option('-i', '--init-params', dest='init_params', type='string', help="Initial parameters separated by commas (no spaces!) Order is: Asini (in lt-s), Orbital period (in days), Pulsar period (in s), T0 (in MJD), Eccentricity, Logitude of periastron (in rad). Example: '5,6.25,0.004909,54905,0.67,-1.51' (NOTE: excluding last two params will result in fitting for a circular orbit.)")
    parser.add_option('-p', '--predict', dest='predict_mjds', type='string', help="MJDs separated by commas at which to predict the spin period (no spaces!)")
    options, args = parser.parse_args()
    
    if options.verbose:
        options.debuglevel = 1
    if options.init_params:
        init_params = np.array([float(x) for x in options.init_params.split(',')])
    else:
        sys.stderr.write("-i, --init-guess option MUST be provided! Exiting!\n")
        sys.exit(1)
    if options.predict_mjds:
        predict_mjds = np.array([float(x) for x in options.predict_mjds.split(',')])
    else:
        predict_mjds = []
        
    main()
