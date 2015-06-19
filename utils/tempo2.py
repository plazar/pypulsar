"""
Some utility functions for making TEMPO2 system calls and parsing output.

Patrick Lazarus, Mar 7, 2014
"""
import tempfile
import subprocess
import os

import numpy as np


def get_resids(parfn, timfn, extra_lines=[], binary=False):
    if extra_lines:
        tmpparhandle, tmpparfn = tempfile.mkstemp(text=True)
        os.close(tmpparhandle)
        tmpparfile = open(tmpparfn, 'w')
        tmpparfile.write(open(parfn).read())
        for ind, offset in dmassplanets:
            tmpparfile.write("\n".join((extra_lines)))
        tmpparfile.close()
    else:
        tmpparfn = parfn
    fmt = r'{bat};;{pre};;{err}'
    if binary:
        fmt += r';;{binphase}'
    pipe = subprocess.Popen(['tempo2', '-output', 'general2', '-f', tmpparfn, \
                                timfn, '-s', fmt+';;\n'], \
                                stdout=subprocess.PIPE)
    stdout, stderr = pipe.communicate()
    datastr = stdout.split("Starting general2 plugin")[1].split(";;\nFinished general2 plugin")[0]
    data = np.fromstring(datastr, sep=';;')
    if binary:
        ncol = 4
    else:
        ncol = 3
    nrow = data.size/ncol
    data.shape = (nrow, ncol) 
    if extra_lines:
        os.remove(tmpparfn)
    return data.T
