#!/usr/bin/env python

import sys
import subprocess
import os
import os.path

debug = False

def get_period_phase_at_mjd(psr, mjd):
    """
    return period and phase of 'psr' at 'mjd'
    
    (Behaviour might be unpredictable if psr cannot
    be found by psrcat.)
    """
    parfilename = "%s.par" % psr
    psrcat_cmd = "psrcat -e %s > %s" % (psr, parfilename)
    if debug: print psrcat_cmd
    psrcat = subprocess.Popen(psrcat_cmd, shell=True)
    psrcat.wait()

    parfile = open(parfilename, 'r')
    lines = [line.split()[-1] for line in parfile.readlines() if line.startswith('PSR')]
    if len(lines): psr = lines[0]
    parfile.close()

    # Check if computer is 32/64-bits
    if os.uname()[4] == 'i686':
        # 32-bit
        bin_dir = "/homes/borg/champion/bin/sigproc-3.5/"
    else:
        # 64-bit
        bin_dir = "/homes/borg/champion/bin/sigproc-3.5-64b/"

    polyco_cmd = "%s %s -mjd %r -par %s > /dev/null" % \
                    (os.path.join(bin_dir, "polyco"), psr, mjd, parfilename)
    if debug: print polyco_cmd
    polyco = subprocess.Popen(polyco_cmd, shell=True)
    polyco.wait()
    
    polyco2period_cmd = "%s %r -p polyco.dat" % \
                            (os.path.join(bin_dir, "polyco2period"), mjd)
    if debug: print polyco2period_cmd
    polyco2period = subprocess.Popen(polyco2period_cmd, shell=True, stdout=subprocess.PIPE)
    polyco2period.wait()
    
    p_phase = polyco2period.stdout.readline().split()
    period = float(p_phase[0])
    phase = float(p_phase[-1])
    polyco2period.stdout.close()
    os.remove(parfilename)

    return period, phase

def main():
    psr = sys.argv[1]
    mjd = float(sys.argv[2])
    period, phase = get_period_phase_at_mjd(psr, mjd)
    print 'PSR %s at MJD %r' % (psr, mjd)
    print 'period: %r\nphase: %r' % (period, phase)

if __name__ == '__main__':
    main()
