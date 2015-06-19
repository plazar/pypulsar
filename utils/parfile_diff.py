"""Plot difference between parfiles.

    Patrick Lazarus, Sept 5, 2014
"""
import sys
import os.path

import numpy as np
import matplotlib.pyplot as plt

import mypolycos as polycos
import psr_utils
from pypulsar import utils


TEL_ID = '3'
FCTR = 1400
MAX_HA = 12

COLOURS = ['r', 'b', 'm', 'c']


def main():
    parfn_ref = sys.argv[1]
    parfns = sys.argv[2:]

    mjds = np.linspace(47000, 48000, 200).astype('f16')
    rots_ref = np.empty(len(mjds))
    rots = np.empty((len(mjds), len(parfns)))

    for ii, mjd in utils.show_progress(enumerate(mjds), width=50, tot=len(mjds)):
        pcos_ref = polycos.create_polycos(parfn_ref, TEL_ID, FCTR, 
                                     np.floor(mjd-1), np.ceil(mjd+1), \
                                     MAX_HA)
        rot = pcos_ref.get_rotation(int(mjd), mjd % 1)
        freq = pcos_ref.get_freq(int(mjd), mjd % 1)
        rots_ref[ii] = np.floor(rot) # Integer rotation number
        # Shift MJD to correspond to integer rotation
        mjd -= (rot % 1)/freq/psr_utils.SECPERDAY
        mjds[ii] = mjd

        for jj, parfn in enumerate(parfns):
            pcos = polycos.create_polycos(parfn, TEL_ID, FCTR, 
                                     np.floor(mjd-1), np.ceil(mjd+1), \
                                     MAX_HA)
            rots[ii,jj] = pcos.get_rotation(int(mjd), mjd % 1)
    fig = plt.figure()
    plt.axhline(0, ls='--', c='k', label=os.path.basename(parfn_ref))
    for jj, parfn in enumerate(parfns):
        c = COLOURS[jj % len(parfns)]
        plt.plot(mjds, rots[:,jj]-rots_ref, c=c, ls='-', lw=2,
                 label=os.path.basename(parfn))
    plt.xlabel("Time (MJD)")
    plt.ylabel("Residuals (revolutions)")
    plt.xlim(mjds.min(), mjds.max())
    plt.legend(loc='best')
    plt.show()


if __name__ == '__main__':
    main()
