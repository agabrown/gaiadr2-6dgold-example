"""
Provide functions that figure out from the astropy.coordinate classes what the rotation matrices between
the various coordinate systems are. This is needed in the context of uncertainty propagation in manner
consistent with the Astropy internals.

Anthony Brown Mar 2019 - Mar 2019
"""

import numpy as np
from astropy.coordinates import Galactic, ICRS
import astropy.units as u

def astropy_icrstogal_rotation_matrix():
    """
    Figure the ICRS to Galatic coordinates rotation matrix.

    Returns
    -------

    The 3x3 matrix representing the rotation from ICRS Cartesian to Galactic Cartesian coordinates.
    """

    # ICRS basis vectors
    basis_icrs = ICRS(ra = [0,np.pi/2,0]*u.rad, dec = [0,0,np.pi/2]*u.rad)

    # Rotate to Galactic
    basis_gal = basis_icrs.transform_to(Galactic())

    # Calculate and return matrix components
    return np.array([np.cos(basis_gal.l.to(u.rad))*np.cos(basis_gal.b.to(u.rad)),
        np.sin(basis_gal.l.to(u.rad))*np.cos(basis_gal.b.to(u.rad)), np.sin(basis_gal.b.to(u.rad))]) 

def astropy_icrstogal_jacobian(icrs_coords, gal_coords):
    """
    For sources with given (ra, dec) and (l,b) calculate the Jacobian for the transformation of the
    astrometric parameters from ICRS to Galactic.

    Parameters
    ----------

    icrs_coords : astropy.coordinates.ICRS
        An astropy.coordinates.ICRS instance containg the list of (ra,dec) coordinates to process.
    gal_coords : astropy.coordinates.Galactic
        An astropy.coordinates.Galactic instance containing the (l,b) coordinates corresponding to the
        ICRS object.

    Returns
    -------

    Array of shape (icrs_coords.ra.size, 5, 5) containing the Jacobians.
    """

    ra = icrs_coords.ra.to(u.rad)
    dec = icrs_coords.dec.to(u.rad)
    l = gal_coords.l.to(u.rad)
    b = gal_coords.b.to(u.rad)
    A = astropy_icrstogal_rotation_matrix()

    picrs = np.array([-np.sin(ra), np.cos(ra), np.zeros_like(ra)])
    qicrs = np.array([-np.cos(ra)*np.sin(dec), -np.sin(ra)*np.sin(dec), np.cos(dec)]) 
    pgal = np.array([-np.sin(l), np.cos(l), np.zeros_like(l)])
    qgal = np.array([-np.cos(l)*np.sin(b), -np.sin(ra)*np.sin(b), np.cos(b)]) 

    J = np.zeros((ra.size, 5, 5))
    J[:,2,2] = 1.0
    for i in range(ra.size):
        G = np.matmul(A, np.stack([picrs[:,i], qicrs[:,i]]).T)
        G = np.matmul(np.stack([pgal[:,i], qgal[:,i]]), G)
        J[i,0:2,0:2] = G
        J[i,3:5,3:5] = G

    return J

if __name__ in ('__main__'): 
    print(astropy_icrstogal_rotation_matrix())

    icrscoo = ICRS(ra = [0,np.pi/2,0,np.pi/3]*u.rad, dec = [0,0,np.pi/2,-np.pi/4]*u.rad)
    galcoo = icrscoo.transform_to(Galactic())
    astropy_icrstogal_jacobian(icrscoo, galcoo)
