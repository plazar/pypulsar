"""
Coordinate conversion and formatting
"""

import psr_utils, slalib

def decstr_to_rad(decstr):
    """
    decstr_to_rad(decstr)

    Convert declination string of the form DDMMSS.SSSS to
    declination in radians
    """
    decstr = str(decstr) # Ensure decstr is string, not int or float    
    sign, d, m, s = parse_decstr(decstr)

    return sign_to_int(sign)*psr_utils.dms_to_rad(float(d), float(m), float(s))

def sign_to_int(sign):
    """
    sign_to_int(sign)

    Convert sign string ('+' or '-') to an integer (1 or -1).
    """

    if sign == '+':
	isign = 1
    elif sign == '-':
	isign = -1
    else:
	raise "ERROR: sign is not '+' or '-' in function sign_to_int."
    return isign
    
def parse_decstr(decstr):
    """
    parse_decstr(decstr)

    Parse declination string of form DDMMSS.SSSS
   
    Returns (sign, d, m, s) where
	sign is '+' or '-'
	d, m, s are degrees, arcmins, arcsecs (strings) respectively
    """
    decstr = str(decstr) # Ensure decstr is string, not int or float    
   
    decl = float(decstr)
    if decl == 0:
	return ('+', '00', '00', '00')
    decl_sign = decl/abs(decl)
    decl = str(abs(decl))
    if '.' in decl:
	decl_split = decl.split('.')
	decl_split[1] = '.%s' % decl_split[1]
    else:
	decl_split = [decl, '']

    decl_pad = decl_split[0].zfill(6)

    if decl_sign == 1:
	sign = '+'
    elif decl_sign == -1:
	sign = '-'
    else:
	raise "ERROR: decl_sign is not +1 or -1 in function parse_decstr"
    d = decl_pad[0:2]
    m = decl_pad[2:4]
    s = '%s%s' % (decl_pad[4:6], decl_split[1])
    return (sign, d, m, s)

def decstr_to_fmdecstr(decstr):
    """
    decstr_to_fmdecstr(decstr)

    Format declination string of form DDMMSS.SSSS
    
    Returns declination string of form +/-DD:MM:SS.SSSS
    """
    decstr = str(decstr) # Ensure decstr is string, not int or float    
    
    return '%s%s:%s:%s' % parse_decstr(decstr)

def decstr_to_deg(decstr):
    """
    decstr_to_deg(decstr)

    Convert declination string of form DDMMSS.SSSS to degrees.
    """
    decstr = str(decstr) # Ensure decstr is string, not int or float    

    return decstr_to_rad(decstr)*psr_utils.RADTODEG

def fmdecstr_to_decstr(fmdecstr):
    """
    fmdecstr_to_decstr(fmdecstr)

    Format declination string of form +/-DD:MM:SS.SSSS
    
    Returns declination string of form DDMMSS.SSSS
    """

    fmdecstr_nocols = fmdecstr.replace(':', '')
    if (fmdecstr_nocols[0] == '-') or (fmdecstr_nocols[0] == '+'):
	sign = fmdecstr_nocols[0]
	fmdecstr_nocols = fmdecstr_nocols[1:]
    else:
	sign = ''
    
    if '.' in fmdecstr_nocols:
	decstr = float(fmdecstr_nocols)
    else:
	decstr = int(fmdecstr_nocols)
    
    return '%s%s' % (sign, decstr)

def rastr_to_rad(rastr):
    """
    rastr_to_rad(rastr)

    Convert right ascension string of form HHMMSS.SSSS to
    radians
    """
    rastr = str(rastr) # Ensure rastr is string, not int or float    
    h, m, s = parse_rastr(rastr)
    return psr_utils.hms_to_rad(float(h), float(m), float(s))

def rastr_to_deg(rastr):
    """
    rastr_to_deg(rastr)

    Convert right ascension string of form HHMMSS.SSSS to
    degrees
    """
    rastr = str(rastr) # Ensure rastr is string, not int or float    
    return rastr_to_rad(rastr)*psr_utils.RADTODEG

def rastr_to_fmrastr(rastr):
    """
    rastr_to_fmrastr(rastr)

    Format right ascension string of from HHMMSS.SSSS.

    Returns right ascension string of form HH:MM:SS.SSSS
    """
    rastr = str(rastr) # Ensure rastr is string, not int or float    
    
    return '%s:%s:%s' % parse_rastr(rastr)

def fmrastr_to_rastr(fmrastr):
    """
    fmrastr_to_rastr(rastr)

    Format right ascension string of form HH:MM:SS.SSSS

    Returns right ascension string of form HHMMSS.SSSS
    """
    
    
    fmrastr_nocols = fmrastr.replace(':', '')
    if '.' in fmrastr_nocols:
	rastr = float(fmrastr_nocols)
    else:
	rastr = int(fmrastr_nocols)
    return '%s' % rastr

def parse_rastr(rastr):
    """
    parse_rastr(rastr)

    Parse right ascension string of form HHMMSS.SSSS

    Returns (h, m, s) where
	h, m, s are hours, minutes, seconds (strings) respectively
    """
    rastr = str(rastr) # Ensure rastr is string, not int or float    

    if float(rastr) == 0:
	return ('00', '00', '00')

    if rastr[0] == '+':
	rastr = rastr[1:]

    if '.' in rastr:
	ra_split = rastr.split('.')
	ra_split[1] = '.%s' % ra_split[1]
    else:
	ra_split = [rastr, '']
    ra_pad = ra_split[0].zfill(6)
    
    h = ra_pad[0:2]
    m = ra_pad[2:4]
    s = '%s%s' % (ra_pad[4:6], ra_split[1])

    return (h, m, s)

def eqdeg_to_galdeg(ra, decl):
    """
    eqdeg_to_galdeg(ra, decl)
    
    Convert right ascension and declination (J2000)in 
    degrees to galactic longitude and latitude in degrees.
    """

    l, b = slalib.sla_eqgal(ra*psr_utils.DEGTORAD, decl*psr_utils.DEGTORAD)

    return (l*psr_utils.RADTODEG, b*psr_utils.RADTODEG)
