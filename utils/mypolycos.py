"""
Taken from Scott Ransom's PRESTO.

Modified to include a function to create polycos from a parfile.
Patrick Lazarus, June 28th, 2009.
"""

import sys
import os
import subprocess
import types
import numpy as Num
import parfile
import psr_utils
import telescopes
import infodata

# Constants
NUMCOEFFS_DEFAULT = 12
SPAN_DEFAULT = 60 # span of each poylco in minutes

class polyco:
    def __init__(self, fileptr):
        line = fileptr.readline()
        if (line==""):
            self.psr = None
        else:
            sl = line.split()
            self.psr = sl[0]
            self.date = sl[1]
            self.UTC = sl[2]
            self.TMIDi = float(sl[3].split(".")[0])
            self.TMIDf = float("0."+sl[3].split(".")[1])
            self.TMID = self.TMIDi+self.TMIDf
            self.DM = float(sl[4])
            if (len(sl)==7):
                self.doppler = float(sl[5])*1e-4
                self.log10rms = float(sl[6])
            else:
                self.log10rms = "-"+sl[-1].split("-")[-1]
                self.doppler = float(sl[-1][:sl[-1].find(self.log10rms)])*1e-4
                self.log10rms = float(self.log10rms)
            sl = fileptr.readline().split()
            self.RPHASE = float(sl[0])
            self.F0 = float(sl[1])
            self.obs = sl[2]
            self.dataspan = int(sl[3])
            self.numcoeff = int(sl[4])
            self.obsfreq = float(sl[5])
            if (len(sl)==7):
                self.binphase = float(sl[6])
            self.coeffs = Num.zeros(self.numcoeff, 'd')
            for linenum in range(self.numcoeff/3):
                sl = fileptr.readline().split()
                self.coeffs[linenum*3+0] = float(sl[0].replace('D', 'E'))
                self.coeffs[linenum*3+1] = float(sl[1].replace('D', 'E'))
                self.coeffs[linenum*3+2] = float(sl[2].replace('D', 'E'))
    
    def phase(self, mjdi, mjdf):
        """
        self.phase(mjdi, mjdf):
            Return the predicted pulsar phase at a given integer and frational MJD.
        """
        return self.rotation(mjdi, mjdf) % 1

    def rotation(self, mjdi, mjdf):
        """
        self.rotation(mjdi, mjdf):
            Return the predicted pulsar (fractional) rotation at a 
            given integer and fractional MJD.
        """
        DT = ((mjdi-self.TMIDi)+(mjdf-self.TMIDf))*1440.0
        phase = self.coeffs[self.numcoeff-1]
        for ii in range(self.numcoeff-1, 0, -1):
            phase = DT*phase + self.coeffs[ii-1]
        phase += self.RPHASE + DT*60.0*self.F0
        return phase 

    def freq(self, mjdi, mjdf):
        """
        self.freq(mjdi, mjdf):
            Return the predicted pulsar spin frequency at a given integer and frational MJD.
        """
        DT = ((mjdi-self.TMIDi)+(mjdf-self.TMIDf))*1440.0
        psrfreq = 0.0
        for ii in range(self.numcoeff-1, 0, -1):
            psrfreq = DT*psrfreq + ii*self.coeffs[ii]
        return self.F0 + psrfreq/60.0


class polycos:
    def __init__(self, filenm="polyco.dat"):
        self.file = filenm
        self.polycos = []
        self.TMIDs = []
        infile = open(filenm, "r")
        tmppoly = polyco(infile)
        psrname = tmppoly.psr
        while(tmppoly.psr):
	    if (len(self.polycos)):
                if (tmppoly.dataspan != self.dataspan):
                    sys.stderr.write("Data span is changing!\n")
                if (tmppoly.psr != psrname):
                    sys.stderr.write("Multiple PSRs in same polycos file!\n")    
            else:
                self.dataspan = tmppoly.dataspan
            self.polycos.append(tmppoly)
            self.TMIDs.append(tmppoly.TMID)
            tmppoly = polyco(infile)
        print "Read %d polycos\n" % len(self.polycos)
        self.TMIDs = Num.asarray(self.TMIDs)
        infile.close()
        self.validrange = 0.5*self.dataspan/1440.0

    def select_polyco(self, mjdi, mjdf):
        """
        self.select_polyco(mjdi, mjdf):
            Return the polyco number that is valid for the specified time.
        """
        goodpoly = Num.argmin(Num.fabs(self.TMIDs-(mjdi+mjdf)))
        if (Num.fabs(self.TMIDs[goodpoly]-(mjdi+mjdf)) > self.validrange):
            sys.stderr.write("Cannot find a valid polyco at %f!\n" % (mjdi+mjdf))
        return goodpoly

    def get_phase(self, mjdi, mjdf):
        """
        self.get_phase(mjdi, mjdf):
            Return the predicted pulsar phase for the specified time.
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return self.polycos[goodpoly].phase(mjdi, mjdf) 

    def get_rotation(self, mjdi, mjdf):
        """
        self.get_rotation(mjdi, mjdf):
            Return the predicted pulsar (fractional) rotation 
            number for the specified time.
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return self.polycos[goodpoly].rotation(mjdi, mjdf) 

    def get_freq(self, mjdi, mjdf):
        """
        self.get_freq(mjdi, mjdf):
            Return the predicted pulsar spin frquency for the specified time.
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return self.polycos[goodpoly].freq(mjdi, mjdf) 

    def get_phs_and_freq(self, mjdi, mjdf):
        """
        self.get_phs_and_freq(mjdi, mjdf):
            Return the predicted pulsar phase and spin frquency for the specified time.
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return (self.polycos[goodpoly].phase(mjdi, mjdf), 
                self.polycos[goodpoly].freq(mjdi, mjdf))
    
    def get_voverc(self, mjdi, mjdf):
        """
        self.get_voverc(mjdi, mjdf):
            Return the (approximate) topocentric v/c for the specified time.
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return self.polycos[goodpoly].doppler


def create_polycos_from_inf(par, infdata):
    """A convenience function to create polycos for the observation
        with info in the given *.inf file.

        Inputs:
            par: parfile's filename, or a parfile object.
            infdata: inffile's filename, or an infodata object.

        Ouput:
            new_polycos: polycos object
    """
    if type(infdata)==types.StringType:
        # assume infdata is a filename
        infdata = infodata.infodata(infdata)
    else:
        # assume infdata is already an infodata.infodata object
        pass
    obslength = (infdata.dt*infdata.N) / psr_utils.SECPERDAY
    telescope_id = telescopes.telescope_to_id[infdata.telescope]
    if telescope_id != 'o' and telescope_id !='@':
        center_freq = infdata.lofreq + (infdata.numchan/2 - 0.5) * \
                                            infdata.chan_width
        if infdata.bary:
            # Data is barycentred, keep max_hour_angle,
            # but set telescope_id to barycentre, '@'
            telescope_id = '@'
    else:
        # optical, X-ray, or gamma-ray data
        center_freq = 0.0

    start_mjd = int(infdata.epoch) 
    end_mjd = int(infdata.epoch+obslength)+1

    return create_polycos(par, telescope_id, center_freq, start_mjd, end_mjd) 


def create_polycos(par, telescope_id, center_freq, start_mjd, end_mjd, \
                    max_hour_angle=None, span=SPAN_DEFAULT, \
                    numcoeffs=NUMCOEFFS_DEFAULT, keep_file=False):
    """Create polycos object from a parfile.
        Inputs:
            par: parfile's filename, or a parfile object.
            telescope_id: The TEMPO 1-character telescope identifier.
            center_freq: The observation's center frequencies in MHz.
            start_mjd: MJD on which the polycos should start.
            end_mjd: MJD until the polycos should extend.
            max_hour_angle: The maximum hour angle as expected by tempo.
                (Default: Use default value chosen for given telescope).
            span: Span of each set of polycos in min.
                (Default: %d min).
            numcoeffs: Number of coefficients to use.
                (Default: %d).
            keep_file: If true do not delete polyco.dat file.
                (Default: delete polyco.dat file).

        Output:
            new_polycos: a polycos object.
    """ % (SPAN_DEFAULT, NUMCOEFFS_DEFAULT)
    if type(par)==types.StringType:
        # assume par is a filename
        par = parfile.psr_par(par)
    else:
        # assume par is already a parfile.psr_par object
        pass
   
    if max_hour_angle is None:
        telescope_name = telescopes.id_to_telescope[telescope_id]
        max_hour_angle = telescopes.telescope_to_maxha[telescope_name]
    
    tzfile = open("tz.in", "w")
    # Default parameters for prediction mode
    tzfile.write("%s %d %d %d %0.5f\n" % (telescope_id, max_hour_angle, \
                            SPAN_DEFAULT, NUMCOEFFS_DEFAULT, center_freq))
    # TEMPO ignores lines 2 and 3 in tz.in file
    tzfile.write("\n\n") 
    tzfile.write("%s %d %d %d %0.5f\n" % (par.PSR, SPAN_DEFAULT, \
                        NUMCOEFFS_DEFAULT, max_hour_angle, center_freq)) 
    tzfile.close()
    tempo = subprocess.Popen("tempo -z -f %s" % par.FILE, shell=True, \
                        stdin=subprocess.PIPE, stdout=subprocess.PIPE, \
                        stderr=subprocess.PIPE)
    (out, err) = tempo.communicate("%d %d\n" % (start_mjd, end_mjd))
    new_polycos = polycos(filenm='polyco.dat')
    # Remove files created by us and by TEMPO
    os.remove("tz.in")
    if not keep_file:
        os.remove("polyco.dat")
    return new_polycos
