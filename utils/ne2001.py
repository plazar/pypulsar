"""
Python interface to NE2001 binary.

Patrick Lazarus, Aug. 16, 2010
"""

import subprocess
import os

if os.uname()[4] == "i686":
    NE2001_LOC = "/homes/borgii/plazar/code/distibutions/NE2001/bin.NE2001"
elif os.uname()[4] == "x86_64":
    NE2001_LOC = "/homes/borgii/plazar/code/distibutions/NE2001_64bit/bin.NE2001"

def get_pulse_broadening(l, b, dm, freq=1):
    """Given galactic longitude (deg), galactic latitude (deg), and
        dispersion measure, return pulse broadening (ms) predicted by
        NE2001 model.
        Optional: frequency in GHz. Scaling uses index of -4.4
    
        NOTE: Only works on computer on which NE2001 was compiled
    """
    # Execute NE2001
    cmd = "./NE2001 %f %f %f 1" % (l, b, dm)
    pipe = subprocess.Popen(cmd, shell=True, cwd=NE2001_LOC, \
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = pipe.communicate()
    for line in output.split('\n'):
        if "PulseBroadening @1GHz" in line:
            broadening = float(line.split()[0])
    # Return scaled broadening
    return broadening*freq**(-4.4)
