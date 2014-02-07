"""
telescopes.py

Various telescope related information.

Patrick Lazarus, July 7th, 2009
"""

# Telescope name to TEMPO observatory code conversion
telescope_to_id = {"GBT": '1', \
                   "Arecibo": '3', \
                   "VLA": '6', \
                   "Parkes": '7', \
                   "Jodrell": '8', \
                   "GB43m": 'a', \
                   "GB 140FT": 'a', \
                   "Nancay": 'f', \
                   "Effelsberg": 'g', \
                   "WSRT": 'i', \
                   "GMRT": 'r', \
                   "Geocenter": 'o', \
                   "Barycenter": '@'}

# TEMPO observatory code to Telescope name conversion
id_to_telescope = {'1': "GBT", \
                   '3': "Arecibo", \
                   '6': "VLA", \
                   '7': "Parkes", \
                   '8': "Jodrell", \
                   'a': "GB43m", \
                   'a': "GB 140FT", \
                   'f': "Nancay", \
                   'g': "Effelsberg", \
                   'i': "WSRT", \
                   'r': "GMRT", \
                   'o': "Geocenter", \
                   '@': "Barycenter"}

# Telescope name to track length (max hour angle) conversion
telescope_to_maxha = {"GBT": 12, \
                   "Arecibo": 3, \
                   "VLA": 6, \
                   "Parkes": 12, \
                   "Jodrell": 12, \
                   "GB43m": 12, \
                   "GB 140FT": 12, \
                   "Nancay": 4, \
                   "Effelsberg": 12, \
                   "WSRT": 12, \
                   "GMRT": 12, \
                   "Geocenter": 12, \
                   "Barycenter": 12}
