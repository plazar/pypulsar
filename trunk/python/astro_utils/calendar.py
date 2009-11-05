"""
calendar.py

A suite of useful calendar/date/day related functions.

Patrick Lazarus, Nov 4th, 2009
"""
import types
import numpy as np


MONTH_TO_NUM = {'January': 1,
                'February': 2,
                'March': 3,
                'April': 4,
                'May': 5,
                'June': 6,
                'July': 7,
                'August': 8,
                'September': 9,
                'October': 10,
                'November': 11,
                'December': 12}

NUM_TO_MONTH = {1: 'January',
                2: 'February',
                3: 'March',
                4: 'April',
                5: 'May',
                6: 'June',
                7: 'July',
                8: 'August',
                9: 'September',
                10: 'October',
                11: 'November',
                12: 'December'}


def JD_to_MJD(JD):
    """Convert Julian Day (JD) to Modified Julian Day (MJD).
    """
    return JD - 2400000.5


def MJD_to_JD(MJD):
    """Convert Modified Julian Day (MJD) to Julian Day (JD).
    """
    return MJD + 2400000.5


def date_to_JD(year, month, day, gregorian=True):
    """Convert calendar date to Julian Day (JD).

        Inputs:
            year: integer
            month:  integer or a string
            day: float
            gregorian:  - True for Gregorian calendar (Default)
                        - False for Julian calendar
    
        (Follow Jean Meeus' Astronomical Algorithms, 2nd Ed., Ch. 7)
    """
    if type(month) == types.StringType:
        month = month_to_num(month)
    
    year = np.atleast_1d(year)
    month = np.atleast_1d(month)
    day = np.atleast_1d(day)
    
    year[month<=2] -= 1
    month[month<=2] += 12
    
    A = np.floor(year/100.0)
    
    if gregorian:
        B = 2 - A + np.floor(A/4.0)
    else:
        B = 0

    JD = np.floor(365.25*(year+4716)) + np.floor(30.6001*(month+1)) + \
            day + B - 1524.5

    if np.any(JD<0.0):
        raise ValueError("This function does not apply for JD < 0.")
    
    return JD.squeeze()


def julian_to_JD(year, month, day):
    """Convert Julian date to Julian Day (JD).

        Inputs:
            year: integer
            month:  integer or a string
            day: float
    
        (Follow Jean Meeus' Astronomical Algorithms, 2nd Ed., Ch. 7)
    """
    return date_to_JD(year, month, day, gregorian=False)


def gregorian_to_JD(year, month, day):
    """Convert Gregorian date to Julian Day (JD).

        Inputs:
            year: integer
            month:  integer or a string
            day: float
    
        (Follow Jean Meeus' Astronomical Algorithms, 2nd Ed., Ch. 7)
    """
    return date_to_JD(year, month, day, gregorian=True)


def month_to_num(month):
    """Return month number given the month name, a string.
    """
    if type(month) != types.StringType:
        raise TypeError("month must be of type string. type(month): %s" % \
                            type(month))
    if month not in MONTH_TO_NUM:
        raise ValueError("Unrecognized month name: %s" % month)
    return MONTH_TO_NUM[month]
