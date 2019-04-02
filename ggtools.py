"""
Provide functions for the building of Gaia DR2 6D Gold sample.

Anthony Brown Mar 2019 - Apr 2019
"""

import numpy as np
import astropy.units as u
from astropy.coordinates import ICRS, Galactic, CartesianDifferential, Galactocentric

# Choice for Sun's phase space coordinates for transformation to Galactocentric reference frame.
_Rsun = 8127.0*u.pc
_zsun = 20.8*u.pc
_vc = 240.0*u.km/u.s
_Usun = 11.1*u.km/u.s
_Vsun = 12.24*u.km/u.s
_Wsun = 7.25*u.km/u.s

def correct_dr2_photometry(Gdr2, GBPdr2, GRPdr2):
    """
    Correct the G, G_BP, G_RP magnitudes as listed in Gaia DR2 according to the recommendations on the
    Gaia DR2 Known Issues pages:
    https://www.cosmos.esa.int/web/gaia/dr2-known-issues#PhotometrySystematicEffectsAndResponseCurves

    Parameters
    ----------

    Gdr2 : float array
        List of G-band magnitudes to process
    GBPdr2 : float array
        List of G_BP-band magnitudes to process
    GRPdr2 : float array
        List of G_RP-band magnitudes to process

    Returns
    -------

    Corrected magnitudes as Gcorr, GBPcorr, GRPcorr
    """

    conditionsG = [(Gdr2 > 6.0) & (Gdr2 <=16), (Gdr2 > 16), (Gdr2 <= 6.0)]
    optionsG = [Gdr2 - 0.0032*(Gdr2-6), Gdr2 - 0.032, -0.047344 + 1.16405*Gdr2 -0.046799*Gdr2**2 +
            0.0035015*Gdr2**3]
    Gcorr = np.select(conditionsG, optionsG)
    GBPcorr = np.where((Gdr2> 2) & (Gdr2< 4), GBPdr2-2.0384+0.95282*Gdr2-0.11018*Gdr2**2, GBPdr2)
    GRPcorr = np.where((Gdr2> 2) & (Gdr2< 3.5), -13.946+14.239*GRPdr2-4.23*GRPdr2**2+0.4532*GRPdr2**3, GRPdr2)

    return Gcorr, GBPcorr, GRPcorr

def transform_to_galactic(icrs_coords):
    """
    For the input astrometry plus radial velocity in the ICRS system calculate the barycentric Galactic
    coordinates as well as Galactocentric coordinates.

    Parameters
    ----------

    icrs_coords: astropy.coordinates.ICRS
        ICRS instance constructed from the 5-parameter Gaia DR2 astrometry and the radial velocity.

    Returns
    -------

    Galactic and Galactocentric objects containing the astrometry in Galactic coordinates, the
    galactocentric Cartesian coordinates, and the galactocentric cylindrical coordinates.
    """

    galactic_coords = icrs_coords.transform_to(Galactic())
    sun_motion = CartesianDifferential(_Usun, _vc+_Vsun, _Wsun)
    galactocentric_cartesian = icrs_coords.transform_to(Galactocentric(galcen_distance=_Rsun, z_sun=_zsun, galcen_v_sun=sun_motion))
    galactocentric_cartesian.set_representation_cls(base='cartesian')
    galactocentric_cylindrical = icrs_coords.transform_to(Galactocentric(galcen_distance=_Rsun, z_sun=_zsun, galcen_v_sun=sun_motion))
    galactocentric_cylindrical.set_representation_cls(base='cylindrical')

    return galactic_coords, galactocentric_cartesian, galactocentric_cylindrical
