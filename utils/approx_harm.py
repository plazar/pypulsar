#!/usr/bin/env python

"""
Approximate the quotient of two periods as a fraction and
an error. This is useful for easily identifying harmonic
relationships.

Patrick Lazarus, June 4, 2011
"""

import sys
import numpy as np

def approx_harm(a, b, maxsteps=20):
    q = [np.nan, np.nan]
    m = [0, 1]
    n = [1, 0]
    x = a
    y = b
    for k in range(2, maxsteps+2):
        if y == 0: 
            break
        q.append(int(x/y))
        tmp = y
        y = x % y
        x = tmp
        m.append(q[k]*m[k-1]+m[k-2])
        n.append(q[k]*n[k-1]+n[k-2])
        #print "%d/%d" % (m[k], n[k])
        if n[k]:
            newfrac = float(m[k])/float(n[k])
            origfrac = float(a)/float(b)
            #print origfrac, newfrac, origfrac-newfrac
            if abs(origfrac-newfrac)<0.01:
                return m[k], n[k]

def output_harm(a, b):
    m, k = approx_harm(a,b)
    origfrac = (float(a)/float(b))

    if m > 9 and k > 9:
        return "%f" % origfrac
    else:
        frac = "%d/%d" % (m, k)
        err = (origfrac - float(m)/float(k))
        if err > 0:
            return "%s + %.2g" % (frac, abs(err))
        elif err < 0:
            return "%s - %.2g" % (frac, abs(err))
        else:
            return "%s" % frac

def main():
    num = float(sys.argv[1])
    den = float(sys.argv[2])
    print output_harm(num, den)

if __name__ == '__main__':
    main()
