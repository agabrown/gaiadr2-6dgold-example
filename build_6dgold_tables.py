"""
Build a "Gaia Gold" catalogue for the 6D sample. Start from a sample downloaded with the
following query:

    select source_id, ra, ra_error, dec, dec_error, parallax, parallax_error, parallax_over_error, pmra,
    pmra_error, pmdec, pmdec_error, ra_dec_corr, ra_parallax_corr, ra_pmra_corr, ra_pmdec_corr,
    dec_parallax_corr, dec_pmra_corr, dec_pmdec_corr, parallax_pmra_corr, parallax_pmdec_corr, 
    pmra_pmdec_corr, radial_velocity, radial_velocity_error,
    phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, bp_rp, g_rp, bp_g, 
    2.5/log(10)*phot_g_mean_flux_over_error as phot_g_mean_mag_error,
    2.5/log(10)*phot_bp_mean_flux_over_error as phot_bp_mean_mag_error,
    2.5/log(10)*phot_rp_mean_flux_over_error as phot_rp_mean_mag_error,
    sqrt(astrometric_chi2_al/(astrometric_n_good_obs_al-5)) as uwe
    from gaiadr2.gaia_source
    where parallax_over_error>5
    and radial_velocity is not null
    and astrometric_params_solved=31
    and rv_nb_transits > 3

For details on the process of building the sample see the Python notebook "Build6DGoldSample.ipynb".

Anthony Brown Mar 2019 - Apr 2019
"""

import numpy as np
import argparse
import os

from astropy.table import Table
from astropy.coordinates import Galactocentric, ICRS, CartesianDifferential
import astropy.units as u

from ruwetools import U0Interpolator
from astropyhacks import *
from ggtools import *

_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    return os.path.join(_ROOT, 'data', path)

def build_catalogue(args):
    """
    Build the catalogue from the input file.

    Parameters
    ----------

    args : list
        Command line arguments

    Returns
    -------

    Nothing
    """

    infile = args['inputFile']

    boubertlist = Table.read("./data/rvscontamination.csv", format="csv")

    poskindata = Table.read(get_data(infile), format='fits')
    sample_size = poskindata['source_id'].size
    print("Number of sources in input table: {0}".format(sample_size))

    rwi = U0Interpolator()
    if (poskindata['phot_g_mean_mag'].size!=sample_size or poskindata['bp_rp'].size!=sample_size):
        print("WARNING: table contains entries without G or (BP-RP).")
    poskindata['ruwe'] = poskindata['uwe']/rwi.get_u0(poskindata['phot_g_mean_mag'], poskindata['bp_rp'])

    poskindata['phot_g_mean_mag_corr'], poskindata['phot_bp_mean_mag_corr'], poskindata['phot_rp_mean_mag_corr'] = \
            correct_dr2_photometry(poskindata['phot_g_mean_mag'], poskindata['phot_bp_mean_mag'], \
            poskindata['phot_rp_mean_mag'])

    lowruweposkin = poskindata['ruwe']<1.40
    reliable_vradposkin = ~np.isin(poskindata['source_id'],boubertlist['source_id'])
    clean_poskin = lowruweposkin & reliable_vradposkin
    print("Number of sources at RUWE<1.4: {0}".format(poskindata['source_id'][lowruweposkin].size))
    print("Number of sources at RUWE<1.4 and with potentially spurious RVs removed: {0}".format(poskindata['source_id'][clean_poskin].size))
    
    icrs_coords = ICRS(ra = poskindata['ra'].to(u.rad),
            dec = poskindata['dec'].to(u.rad),
            distance = (1000/poskindata['parallax'])*u.pc/u.mas,
            pm_ra_cosdec = poskindata['pmra'],
            pm_dec = poskindata['pmdec'],
            radial_velocity = poskindata['radial_velocity'])

    galactic_coords, galactocentric_cartesian, galactocentric_cylindrical = transform_to_galactic(icrs_coords)

    poskindata['l'] = galactic_coords.l.to(u.deg)
    poskindata['b'] = galactic_coords.b.to(u.deg)
    poskindata['pml'] = galactic_coords.pm_l_cosb
    poskindata['pmb'] = galactic_coords.pm_b

    poskindata['x_gc'] = galactocentric_cartesian.x
    poskindata['y_gc'] = galactocentric_cartesian.y
    poskindata['z_gc'] = galactocentric_cartesian.z
    poskindata['v_x_gc'] = galactocentric_cartesian.v_x
    poskindata['v_y_gc'] = galactocentric_cartesian.v_y
    poskindata['v_z_gc'] = galactocentric_cartesian.v_z

    #Convert Cylindrical into conventional units (km/s for the velocities, making v_phi positive along
    #the direction of Galactic rotation).
    poskindata['R_gc'] = galactocentric_cylindrical.rho
    phi = galactocentric_cylindrical.phi.to(u.deg)
    poskindata['phi_gc'] = np.where(phi<0*u.deg, phi+360*u.deg, phi.to(u.deg))*u.deg
    poskindata['v_R_gc'] = galactocentric_cylindrical.d_rho.to(u.km/u.s)
    # In the literature vphi is calculated for a left-handed coordinate system! 
    # This is for the convenience of having postive values of vphi at the position of the sun.
    poskindata['v_phi_gc'] = -(galactocentric_cylindrical.d_phi.to(u.rad/u.yr)/u.rad * galactocentric_cylindrical.rho).to(u.km/u.s)

    poskindata[clean_poskin].write("GaiaDR2-6DGold.fits", format="fits")

def parseCommandLineArguments():
    """
    Set up command line parsing.
    """
    parser = argparse.ArgumentParser("Build the Gaia Gold 6D catalogue.")
    parser.add_argument('inputFile', type=str, help="""FITS file from which to build the catalogue""")
    args = vars(parser.parse_args())
    return args

if __name__ in ('__main__'):
    args=parseCommandLineArguments()
    build_catalogue(args)
