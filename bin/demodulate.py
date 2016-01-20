import argparse
import tempfile
import os.path
import types
import sys

import numpy as np

import polycos
import infodata
import psr_utils

from pypulsar.formats import datfile
from pypulsar import utils


def create_polycos_from_inf(par, inf):
    """A convenience function to create polycos for the observation
        with info in the given *.inf file.

        Inputs:
            par: parfile's filename, or a parfile object.
            inf: inffile's filename, or an infodata object.

        Ouput:
            new_polycos: polycos object
    """
    if type(inf) == types.StringType:
        # assume inf is a filename
        inf = infodata.infodata(inf)
    else:
        # assume inf is already an infodata.infodata object
        pass
    obslength = (inf.dt*inf.N) / psr_utils.SECPERDAY
    telescope_id = polycos.telescope_to_id[inf.telescope]
    if telescope_id != 'o' and telescope_id != '@':
        center_freq = inf.lofreq + (inf.numchan/2 - 0.5) * inf.chan_width
        if inf.bary:
            # Data is barycentred, keep max_hour_angle,
            # but set telescope_id to barycentre, '@'
            telescope_id = '@'
    else:
        # optical, X-ray, or gamma-ray data
        center_freq = 0.0

    start_mjd = int(inf.epoch) 
    end_mjd = int(inf.epoch+obslength)+1

    return polycos.create_polycos(par, telescope_id, center_freq, start_mjd, end_mjd) 


def create_parfile(inparfn, inf):
    outfd, outfn = tempfile.mkstemp(suffix='.par', dir=os.getcwd(), text=True)
    outff = os.fdopen(outfd, 'w')
    outff.write("RAJ %s\n" % inf.RA)
    outff.write("DECJ %s\n" % inf.DEC)
    # Set folding frequency to be 1000th of sampling rate (i.e. 1/dt)
    # for it to be comparable with pulsar periods (TEMPO polycos format
    # doesn't allow for enough digits for ~50us periods)
    outff.write("F0 %.15f\n" % (0.001/inf.dt))
    outff.write("F1 0\n")
    outff.write("DM 0\n")
    outff.write("PEPOCH %.15f\n" % inf.epoch)
    outff.write("POSEPOCH %.15f\n" % inf.epoch)
    outff.write("TZRMJD %.15f\n" % inf.epoch)
    outff.write("TZRSITE @\n")  # Barycentre
    outff.write("TZRFREQ %.5f\n" % (inf.lofreq+0.5*inf.BW))

    with open(inparfn, 'r') as inff:
        for line in inff:
            line = line.strip()
            split = line.split()
            key = split[0]
            if key not in ["F", "F0", "F1", "F2", "F3", "F4", "F5", "F6",
                           "P", "P0", "P1", "P2", "P3", "P4", "P5", "P6",
                           "RAJ", "DECJ", "ELAT", "ELONG", "LAMBDA", "BETA",
                           "LAMBDA", "BETA", "RA_RAD", "DEC_RAD", "PMRA", 
                           "PMDEC", "PEPOCH", "POSEPOCH"]:
                outff.write(" ".join(split[0:2])+'\n')
    outff.close()
    return outfn


def main():
    indat = datfile.Datfile(args.datfile)
    data = indat.read_all()
    # Construct parfile to use based on input parfile
    # (folding frequency is 1/1000 of sampling time)
    parfn = create_parfile(args.parfile, indat.inf)
    # Create polycos
    pcos = create_polycos_from_inf(parfn, indat.inf) 

    outdata = []
    dphase = 0.0
    imjd = np.floor(indat.inf.epoch)
    fmjd = indat.inf.epoch % 1

    rot0 = pcos.get_rotation(imjd, fmjd)
    #print "MJD at file start: %.15f" % indat.inf.epoch
    #print "imjd: %05d; fmjd: %.15f" % (imjd, fmjd)
    #print "Rotation at file start: %.15f" % rot0
    samp_in_day = indat.inf.dt/psr_utils.SECPERDAY

    nremoved = 0
    nadded = 0
    last = False
    mjdend = indat.inf.epoch+indat.inf.dt*indat.inf.N/psr_utils.SECPERDAY
    status = 0
    totsamps = indat.inf.N
    startsamp = 0

    while not last:
        igoodpoly = pcos.select_polyco(imjd, fmjd)
        pco = pcos.polycos[igoodpoly]
        pcoend = pco.TMID + pcos.validrange
        nsamp = min(indat.inf.N, int((pcoend - (imjd+fmjd))*psr_utils.SECPERDAY/indat.inf.dt))
        last = (pcoend > mjdend)
        for idatsamp in xrange(startsamp, startsamp+nsamp):
            # Update integer and fractional part of MJD
            # accounting for crossing over to a new day
            newday = (fmjd > 1.0)
            imjd += newday
            fmjd -= newday

            # Compute the sample number in the pulsar's frame using the polycos
            new_rot = pco.rotation(imjd, fmjd)
            ipsrsamp = (new_rot-rot0)*1000.0  # multiply by 1000 because period
                                              # is 1000x sample time
            # Calculate difference in sample number in observation
            # and pulsar frames
            dsamp = (idatsamp+dphase)-ipsrsamp
            #print "Data samp#: %d; dphase: %d; imjd: %05d; fmjd: %.15f; Rot#: %.6f; dRot: %.10f; PSR samp#: %.3f; diff: %.10f" % \
            #      (idatsamp, dphase, imjd, fmjd, new_rot, (new_rot-rot0), ipsrsamp, dsamp)
            if dsamp > 0.5:
                # Drop a sample
                dphase -= 1
                nremoved += 1
            elif dsamp < -0.5:
                # Add a sample
                dphase += 1
                outdata.append(data[idatsamp-1])
                # And now the extra sample is a duplicate
                outdata.append(data[idatsamp-1])
                nadded += 1
            else:
                outdata.append(data[idatsamp-1])
            #print "dphase: %.15f; Nadd: %d; Nremove: %d" % (dphase, nadded, nremoved)
            fmjd += samp_in_day

            # Show status
            newstatus = int(100.0*idatsamp/totsamps)
            if newstatus > status:
                status = newstatus
                sys.stdout.write(" %d %%\r")
                sys.stdout.flush()
        startsamp += nsamp
    sys.stdout.write("\nDone\n")
    print "Number of samples removed: %d" % nremoved
    print "Number of samples added: %d" % nadded
    outdata = np.array(outdata)
    outdata.astype('float32')
    outdata.tofile(args.datfile[:-4]+"_demod.dat")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("datfile",
                        help="PRESTO *.dat file to demodulate. (Corresponding "
                             "*.inf file must be in same directory).")
    parser.add_argument("-f", "--parfile", dest='parfile', required=True,
                        help="Parfile to use to de-modulate orbit.")
    args = parser.parse_args()
    main()
