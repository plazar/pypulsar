"""
sextant.py

A suite of useful coordinate related functions.

Patrick Lazarus, Jan. 31, 2010
"""
import types
import numpy as np

import protractor

def ha_from_lst(lst, ra):
    """Given local sidereal time (lst) and right ascension (ra),
        both in same units, (degrees, radians, or decimal hours)
        calculate the local hour angle.

        Return value is in same units as input arguments.
    """
    hourangle = lst - ra
    return hourangle


def ha_from_mjdlon(mjd, lon, ra):
    """Given Modified Julian Day (mjd), longitude (lon), right ascension (ra)
        where longitude, measured positive west and negative east from
        Greenwich, is given in degrees and right ascension is also given
        in degrees.

        Return value is in degrees.
    """
    mst_deg = clock.MJD_to_mstUT_deg(mjd)
    hourangle = mst_deg - lon - ra
    return hourangle


def equatorial_to_ecliptic(ra, decl, input="sexigesimal", output="deg", \
                                J2000=True):
    """Given right ascension and declination (provided in units of 'input')
        return ecliptic coordinates (longitude, latitude) in units of 'output'.

        If J2000 is true, assume input coords are in J2000 equinox,
            otherwise assume input coords are in B1950 equinox.
    
        Possible values for input and output are "sexigesimal", "deg" and "rad".
    """
    if J2000:
        obliquity = 0.409092804 # radians
    else:
        obliquity = 0.409206212 # radians

    # Convert equatorial coords to radians
    if input == "sexigesimal":
        ra = protractor.convert(ra, "hmsstr", "rad")
        decl = protractor.convert(decl, "dmsstr", "rad")
    else:
        ra = protractor.convert(ra, input, "rad")
        decl = protractor.convert(decl, input, "rad")

    # Do the conversion
    lon = np.arctan2(np.sin(ra)*np.cos(obliquity) + \
                        np.tan(decl)*np.sin(obliquity),
                        np.cos(ra))
    lat = np.arcsin(np.sin(decl)*np.cos(obliquity) - \
                        np.cos(decl)*np.sin(obliquity)*np.sin(ra))
    
    # Ensure radian values are between 0 and 2pi
    lon = np.mod(lon, np.pi*2)
    lat = np.mod(lat, np.pi*2)
    
    if output == "sexigesimal":
        output = "dmsstr"
    lon = protractor.convert(lon, "rad", output)
    lat = protractor.convert(lat, "rad", output)

    return (lon, lat)


def ecliptic_to_equatorial(lon, lat, input="deg", output="sexigesimal", \
                                J2000=True):
    """Given ecliptic longitude and latitude (provided in units of 'input')
        return equatorial coordinates (right ascension, declination)
        in units of 'output'.

        If J2000 is true, assume input coords are in J2000 equinox,
            otherwise assume input coords are in B1950 equinox.
    
        Possible values for input and output are "sexigesimal", "deg" and "rad".
    """
    if J2000:
        obliquity = 0.409092804 # radians
    else:
        obliquity = 0.409206212 # radians

    # Convert ecliptic coords to radians
    if input == "sexigesimal":
        input = "dmsstr"
    lon = protractor.convert(lon, input, "rad")
    lat = protractor.convert(lat, input, "rad")

    # Do the conversion
    ra = np.arctan2(np.sin(lon)*np.cos(obliquity) - \
                    np.tan(lat)*np.sin(obliquity), \
                    np.cos(lon))
    decl = np.arcsin(np.sin(lat)*np.cos(obliquity) + \
                    np.cos(lat)*np.sin(obliquity)*np.sin(lon))
   
    # Ensure radian values are between 0 and 2pi
    ra = np.mod(ra, np.pi*2)
    decl = np.mod(decl, np.pi*2)
    
    if output == "sexigesimal":
        ra = protractor.convert(ra, "rad", "hmsstr")
        decl = protractor.convert(decl, "rad", "dmsstr")
    else:
        ra = protractor.convert(ra, "rad", output)
        decl = protractor.convert(decl, "rad", output)

    return (ra, decl)


def hadec_to_altaz(ha, decl, obslat, input="sexigesimal", output="deg"):
    """Given hour angle, declination (provided in units of 'input') 
        and observer latitude (in degrees) return local horizontal coordinates
        (altitude and azimuth, in units of 'output').
    
        Possible values for input and output are "sexigesimal", "deg" and "rad".
    """
    # Convert equatorial coords to radians
    if input == "sexigesimal":
        ha = protractor.convert(ha, "hmsstr", "rad")
        decl = protractor.convert(decl, "dmsstr", "rad")
    else:
        ha = protractor.convert(ha, input, "rad")
        decl = protractor.convert(decl, input, "rad")

    # Do the conversion
    az = np.arctan2(np.sin(ha), np.cos(ha)*np.sin(obslat) - \
                        np.tan(decl)*np.cos(obslat))
    alt = np.arcsin(np.sin(obslat)*np.sin(decl) + \
                        np.cos(obslat)*np.cos(decl)*np.cos(ha))

    # Ensure radian values are between 0 and 2pi
    az = np.mod(az, np.pi*2)
    alt = np.mod(alt, np.pi*2)

    # Convert output values to desired units
    if output == "sexigesimal":
        output = "dmsstr"
    az = protractor.convert(az, "rad", output)
    alt = protractor.convert(alt, "rad", output)

    return (alt, az)


def altaz_to_hadec(alt, az, obslat, input="deg", output="sexigesimal"):
    """Given altitude, azimuth angle (provided in units of 'input')
        and observer latitude (in degrees) return hour angle and 
        declination (in units of 'output').

        Possible values for input and output are "sexigesimal", "deg" and "rad".
    """
    # Convert input args to radians
    if input == "sexigesimal":
        input = "dmsstr"
    alt = protractor.convert(alt, input, "rad")
    az = protractor.convert(az, input, "rad")
    
    # Do the conversion
    ha = np.arctan2(np.sin(az), np.cos(az)*np.sin(obslat) + \
                    np.tan(alt)*np.cos(obslat))
    decl = np.arcsin(np.sin(obslat)*np.sin(alt) - \
                    np.cos(obslat)*np.cos(alt)*np.cos(az))

    # Ensure radian values are between 0 and 2pi
    ha = np.mod(ha, np.pi*2)
    decl = np.mod(decl, np.pi*2)

    # Convert output values to desired units
    if output == "sexigesimal":
        ha = protractor.convert(ha, "rad", "hmsstr")
        decl = protractor.convert(decl, "rad", "dmsstr")
    else:
        ha = protractor.convert(ha, "rad", output)
        decl = protractor.convert(decl, "rad", output)

    return (ha, decl)


def equatorial_to_galactic(ra, decl, input="sexigesimal", output="deg", \
                            equinox="J2000"):
    """Given right ascension and declination (in units of 'input') convert
        to galactic longitude and latitude (returned in units of 'output').

        Possible values for input and output are "sexigesimal", "deg" and "rad".
        'equinox' is the standard equinox used for the input equatorial
            coordinates. It can be "J2000" (the default), or "B1950".
    """
    # Convert equatorial coords to radians
    if input == "sexigesimal":
        ra = protractor.convert(ra, "hmsstr", "rad")
        decl = protractor.convert(decl, "dmsstr", "rad")
    else:
        ra = protractor.convert(ra, input, "rad")
        decl = protractor.convert(decl, input, "rad")

    # Define galactic north pole
    if equinox == "J2000":
        raise NotImplementedError("J2000 values are incorrect.")
        ra_north = 3.36603297 # radians
        decl_north = 0.47347819 # radians
    elif equinox == "B1950":
        ra_north = 3.35539549 # radians
        decl_north = 0.478220215 # radians
    else:
        raise ValueError("Unrecognized 'equinox': %s" % equinox)

    # Do the conversion
    x = np.arctan2(np.sin(ra_north-ra), np.cos(ra_north-ra)*np.sin(decl_north) - \
                    np.tan(decl)*np.cos(decl_north))
    l = 5.28834763 - x # 303 deg = 5.28834763 rad (origin of galactic coords)
    b = np.arcsin(np.sin(decl)*np.sin(decl_north) + \
                    np.cos(decl)*np.cos(decl_north)*np.cos(ra_north-ra))

    # Ensure radian values are between 0 and 2pi
    l = np.mod(l, np.pi*2)
    b = np.mod(b, np.pi*2)

    # Convert output values to desired units
    if output == "sexigesimal":
        output = "dmsstr"
    l = protractor.convert(l, "rad", output)
    b = protractor.convert(b, "rad", output)

    return (l, b)
