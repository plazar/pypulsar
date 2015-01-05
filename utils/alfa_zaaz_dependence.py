import numpy as np
import matplotlib.pyplot as plt

"""
Zenith angle dependence of ALFA's gain, Tsys, and SEFD.

Data are taken from tarball http://www.naic.edu/alfa/performance/ALFA_POLY_FITS.tar.gz

at http://www.naic.edu/alfa/performance/spectral-scans/wide/

Specifically:
SEFD_Vs_ZA_beam0_olddata_fit.parameters
Gain_Vs_ZA_beam0_olddata_fit.parameters
Tsys_Vs_ZA_beam0_olddata_fit.parameters

Patrick Lazarus, Nov 21, 2014
"""

# Added default value before other params
GAIN_PARAM_STR = """
 10.4
 0 1  5.0  19.3700008  10.043704  10.043704 11 15 41  5.9939723 -0.624729395
  1.52758908 -1.08500731  0.606789947 -1.49469185  0.152855217 -1.87550592
 -0.156861529 -2.22461319 -0.398988336  4.2598381 -0.391409189  0.685782075
  0.792036533 -1.31411183  0.603479087 -0.371651351 -1.30490589  0.889832795
 -0.593093336  0.0949792564  1.83947074 -0.741901636  0.333228111  0.323233545
 -2.47698593  0.539871395  0.283156157 -0.988350868  3.07428741  0.213247508
 -1.73438001  1.72857463 -2.91462374 -2.96988988  4.98494482  2.21380353
 -3.12255979 -0.691958249  0.777421355  0.00988082867 -15.
"""
# Added default value before other params
SEFD_PARAM_STR = """
 3
 0 1  5.0  19.3700008  10.043704  10.043704 11 5 21  2.07651114  0.0696394295
  0.962545931  0.0991852432  0.751455009  0.1668275  0.455828071  0.204119235
 -0.117904358  0.094586201 -0.907949626  1.07005715  0.0577052683 -0.239431992
  0.0185407307  0.186046168  0.127920657 -0.0259651244 -0.203498781
 -0.0168917663  0.0998328701  0.0140674142  7.
"""
# Added default value before other params
TSYS_PARAM_STR = """
 29
 0 1  5.0  19.3700008  10.043704  10.043704 6 2 10  28.4584408  0.627815545
  26.8757477  1.04016066 -15.9114399  1.35548031 -5.35760641  0.422170252
  6.97873116 -0.0233611483  0.176407114  18.
"""


def parse_params(paramstr):
    vals = [float(x) for x in paramstr.strip().split()]
    params = {}
    params['default'] = float(vals[0])
    params['beam_id'] = int(vals[1])
    params['pol_id'] = int(vals[2])
    params['start_za'] = vals[3]
    params['stop_za'] = vals[4]
    params['ref_za'] = vals[5]
    params['halfspan_za'] = vals[6]
    num_poly_coeffs = int(vals[7])
    num_harm_coeffs = int(vals[8])
    num_tot_coeffs = int(vals[9])
    params['poly_coeffs'] = vals[10:10+num_poly_coeffs]
    params['cos_coeffs'] = vals[10+num_poly_coeffs:10+num_tot_coeffs:2]
    params['sin_coeffs'] = vals[1+10+num_poly_coeffs:10+num_tot_coeffs:2]

    return params


def zaaz_func_factory(params):
    """Return function that is dependent on zenith angle and azimuth
    
        Inputs:
            params: A dictionary of parameters with the following keys:
                beam_id - ignored
                pol_id - ignored
                start_za
                stop_za
                ref_za
                halfspan_za
                poly_coeffs
                cos_coeffs
                sin_coeffs

        Output:
            zaaz_func: A function that takes za,az (in degrees) as arguments
                and returns a value. What that value is depeneds on the input
                parameters.
    """
    def zaaz_func(za, az=None):
        scaledza = (za - params['ref_za'])/params['halfspan_za']

        poly = np.polyval(params['poly_coeffs'][::-1], scaledza)
        angles = scaledza[:,np.newaxis]*np.arange(1, 1+len(params['cos_coeffs']))*np.pi/2.0
        cos = params['cos_coeffs']*np.cos(angles)
        sin = params['sin_coeffs']*np.sin(angles)
        vals = poly+np.sum(cos, axis=1)+np.sum(sin, axis=1)
        vals[(za < params['start_za']) | (za > params['stop_za'])] = params['default']
        return vals
    return zaaz_func


gain_params = parse_params(GAIN_PARAM_STR)
gain = zaaz_func_factory(gain_params)

sefd_params = parse_params(SEFD_PARAM_STR)
sefd = zaaz_func_factory(sefd_params)

tsys_params = parse_params(TSYS_PARAM_STR)
tsys = zaaz_func_factory(tsys_params)


def main():
    za = np.linspace(0, gain_params['stop_za'], 1000)

    plt.figure()
    plt.plot(za, gain(za), 'k-', label="Beam 0")
    plt.plot(za, gain(za)/10.4*8.2, 'k:', label="Beams 1-6")
    plt.xlabel("Zenith Angle (deg)")
    plt.ylabel("Gain (K/Jy)")
    plt.legend(loc='best')

    plt.figure()
    plt.plot(za, sefd(za), 'k-')
    plt.plot(za, tsys(za)/gain(za), 'r-')
    plt.xlabel("Zenith Angle (deg)")
    plt.ylabel("SEFD (Jy)")

    plt.figure()
    plt.plot(za, tsys(za), 'k-')
    plt.xlabel("Zenith Angle (deg)")
    plt.ylabel("Tsys (K)")

    plt.show()

if __name__ == '__main__':
    main()
