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

    if args.outname is None:
        outname = indat.basefn
    else:
        outname = args.outname
    if os.path.exists(outname+".dat") and not args.force:
        raise ValueError("Output file (%s) already exists!" % (outname+".dat"))
    if os.path.exists(outname+".inf") and not args.force:
        raise ValueError("Output file (%s) already exists!" % (outname+".inf"))

    # Construct parfile to use based on input parfile
    # (folding frequency is 1/1000 of sampling time)
    parfn = create_parfile(args.parfile, indat.inf)
    # Create polycos
    pcos = create_polycos_from_inf(parfn, indat.inf) 

    dphase = 0.0
    imjd = np.floor(indat.inf.epoch)
    fmjd = np.float128(indat.inf.epoch) % 1

    rot0 = pcos.get_rotation(imjd, fmjd)
    samp_in_day = np.float128(indat.inf.dt)/psr_utils.SECPERDAY

    nremoved = 0
    nadded = 0
    last = False
    mjdend = indat.inf.epoch+indat.inf.dt*indat.inf.N/psr_utils.SECPERDAY
    status = 0
    totsamps = indat.inf.N
    startsamp = 0
    idrop = []
    iadd = []

    idatsamp = 0
    igoodpoly = pcos.select_polyco(imjd, fmjd)
    while not last:
        pco = pcos.polycos[igoodpoly]
        pcoend = pco.TMID + pcos.validrange
        sampstep = 1
        last_add_drop = startsamp
        nsamp = min(indat.inf.N-startsamp,
                    int((pcoend - (imjd+fmjd))*psr_utils.SECPERDAY/indat.inf.dt))
        last = (pcoend > mjdend)
        while idatsamp < (startsamp+nsamp):
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
            if dsamp > 0.5:
                if sampstep == 1:
                    idrop.append(idatsamp)
                    sampstep = int(0.1*(idatsamp - last_add_drop))
                    last_add_drop = idatsamp
                else:
                    steps = np.arange(0, sampstep)
                    fmjds = fmjd - steps*samp_in_day
                    rots = pco.rotation(imjd, fmjds)
                    ipsrsamps = (rots - rot0)*1000.0
                    dsamps = (idatsamp-steps+dphase)-ipsrsamps
                    todrop = np.argmax(dsamps <= 0.5)
                    dropsamp = idatsamp-(todrop-1)
                    idrop.append(dropsamp)
                    sampstep = int(0.1*(dropsamp - last_add_drop))
                    last_add_drop = dropsamp
                # Drop a sample
                dphase -= 1
                nremoved += 1
            elif dsamp < -0.5:
                if sampstep == 1:
                    iadd.append(idatsamp)
                    sampstep = int(0.1*(idatsamp - last_add_drop))
                    last_add_drop = idatsamp
                else:
                    steps = np.arange(0, sampstep)
                    fmjds = fmjd - steps*samp_in_day
                    rots = pco.rotation(imjd, fmjds)
                    ipsrsamps = (rots - rot0)*1000.0
                    dsamps = (idatsamp-steps+dphase)-ipsrsamps
                    toadd = np.argmax(dsamps >= -0.5)
                    addsamp = idatsamp-(toadd-1)
                    iadd.append(addsamp)
                    sampstep = int(0.1*(addsamp - last_add_drop))
                    last_add_drop = addsamp
                # Add a sample
                dphase += 1
                nadded += 1
            else:
                pass
            if sampstep < 1:
                sampstep = 1
            # Prepare for next iteration
            if (idatsamp+sampstep) > (startsamp+nsamp):
                sampstep = (startsamp+nsamp)-idatsamp
            idatsamp += sampstep
            fmjd += samp_in_day * sampstep

            # Show status
            newstatus = int(100.0*idatsamp/totsamps)
            if newstatus > status:
                status = newstatus
                sys.stdout.write(" %d %%\r" % status)
                sys.stdout.flush()
        startsamp += nsamp
        igoodpoly += 1
    sys.stdout.write("\nDone\n")
    print "Number of samples removed: %d" % nremoved
    print "Number of samples added: %d" % nadded

    samps = np.concatenate((idrop, iadd)).astype('int32')
    isdrops = np.zeros_like(samps, dtype='int8')
    isdrops[:len(idrop)] = 1
    isort = np.argsort(samps)
    samps = samps[isort]
    isdrops = isdrops[isort]

    indat.rewind()  # Rewind file, just in case
    with open(outname+".dat", 'w') as outff:
        for ind, isdrop in zip(samps, isdrops):
            data = indat.read_to(ind)
            if isdrop:
                data[:-1].tofile(outff)  # drop last sample
            else:
                data.tofile(outff)
                data[-1:].tofile(outff)  # duplicate last sample
        data = indat.read_to(-1)  # Read rest of file
        data.tofile(outff)
    
    indat.inf.deorbited = True
    indat.inf.to_file(outname+".inf")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("datfile",
                        help="PRESTO *.dat file to demodulate. (Corresponding "
                             "*.inf file must be in same directory).")
    parser.add_argument("-f", "--parfile", dest='parfile', required=True,
                        help="Parfile to use to de-modulate orbit.")
    parser.add_argument("-o", "--outname", dest='outname', default=None,
                        help="Base name (i.e. no extension) of the output "
                             "time series file. (Default: <input base name>_demod)")
    parser.add_argument("--force", dest='force', action='store_true',
                        help="Overwrite any existing files with names of "
                             "output files. (Default: Don't)")
    args = parser.parse_args()
    main()
