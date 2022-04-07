#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Column names from RM-tools to catalogue

@author: cvaneck
"""

# Each column is a tuple:
# output name
# variable type (int/float/string/etc)
# source, can be synth (from RMtools_1d), cat (from input Stokes I catalog),
#     header (from FITS header), meta (from meta data)
# dict key or column name in pipeline
# unit (string)

import astropy.units as u
from astropy.units import cds
cds.enable()

columns = [
    # Relevant quantities from the source finder:
    ('ra', float, 'cat', 'RA', u.deg),
    ('ra_err', float, 'cat', 'E_RA', u.arcsec),
    ('dec', float, 'cat', 'Dec', u.deg),
    ('dec_err', float, 'cat', 'E_Dec', u.arcsec),
    ('total_I_flux', float, 'cat', 'Total_flux_Gaussian', u.mJy),
    ('total_I_flux_err', float, 'cat', 'E_Total_flux_Gaussian', u.mJy),
    ('peak_I_flux', float, 'cat', 'Peak_flux', u.mJy/u.beam),
    ('peak_I_flux_err', float, 'cat', 'E_Peak_flux', u.mJy/u.beam),
    ('maj', float, 'cat', 'Maj', u.arcsec),
    ('maj_err', float, 'cat', 'E_Maj', u.arcsec),
    ('min', float, 'cat', 'Min', u.arcsec),
    ('min_err', float, 'cat', 'E_Min', u.arcsec),
    ('pa', float, 'cat', 'PA', u.deg),
    ('pa_err', float, 'cat', 'E_PA', u.deg),
    ('dc_maj', float, 'cat', 'DC_Maj', u.arcsec),
    ('dc_maj_err', float, 'cat', 'E_DC_Maj', u.arcsec),
    ('dc_min', float, 'cat', 'DC_Min', u.arcsec),
    ('dc_min_err', float, 'cat', 'E_DC_Min', u.arcsec),
    ('dc_pa', float, 'cat', 'DC_PA', u.deg),
    ('dc_pa_err', float, 'cat', 'E_DC_PA', u.deg),
    ('stokesI_err', float, 'cat', 'Noise', u.mJy/u.beam),
    ('beamdist', float, 'cat', 'Separation_Tile_Centre', u.deg),
    ('N_Gaus', int, 'cat', 'N_Gaus', u.dimensionless_unscaled),
    ('cat_id', str, 'cat', 'Gaussian_ID', None),
    ('source_id', str, 'cat', 'Source_ID', None),
    ('tile_id', str, 'cat', 'Tile_ID', None),
    ('sbid', str, 'cat', 'SBID', None),
    ('start_time', float, 'cat', 'Obs_Start_Time', cds.MJD),
    ('separation_tile_centre', float, 'cat', 'Separation_Tile_Centre', u.deg),
    ('s_code', str, 'cat', 'S_Code', None),


    # Important quantities from the RMsynthesis
    ('rm', float, 'synth', 'phiPeakPIfit_rm2', u.rad/u.m**2),
    ('rm_err', float, 'synth', 'dPhiPeakPIfit_rm2', u.rad/u.m**2),
    ('polint', float, 'synth', 'ampPeakPIfitEff', u.Jy/u.beam),
    ('polint_err', float, 'synth', 'dAmpPeakPIchan', u.Jy/u.beam),
    ('stokesQ', float, 'synth', 'peakFDFrealFit', u.Jy/u.beam),
    ('stokesU', float, 'synth', 'peakFDFimagFit', u.Jy/u.beam),
    ('polangle', float, 'synth', 'polAngleFit_deg', u.deg),
    ('polangle_err', float, 'synth', 'dPolAngleFit_deg', u.deg),
    ('derot_polangle', float, 'synth', 'polAngle0Fit_deg', u.deg),
    ('derot_polangle_err', float, 'synth', 'dPolAngle0Fit_deg', u.deg),
    ('fracpol', float, 'synth', 'fracPol', u.dimensionless_unscaled),
    ('reffreq_pol', float, 'synth', 'freq0_Hz', u.Hz),
    ('reffreq_beam', float, 'synth', 'freq0_Hz', u.Hz),
    ('reffreq_I', float, 'synth', 'freq0_Hz', u.Hz),
    ('fdf_noise_th', float, 'synth', 'dFDFth', u.Jy/u.beam),
    ('rmsf_fwhm', float, 'synth', 'fwhmRMSF', u.rad/u.m**2),
    ('refwave_sq_pol', float, 'synth', 'lam0Sq_m2', u.m**2),
    ('stokesI', float, 'synth', 'Ifreq0', u.Jy/u.beam),
    ('stokes_I_fit_flag', int, 'synth', 'IfitStat', u.dimensionless_unscaled),
    ('snr_polint', float, 'synth', 'snrPIfit', u.dimensionless_unscaled),
    ('minfreq', float, 'synth', 'min_freq', u.Hz),
    ('maxfreq', float, 'synth', 'max_freq', u.Hz),
    ('channelwidth', float, 'synth', 'median_channel_width', u.Hz),
    ('Nchan', int, 'synth', 'N_channels', u.dimensionless_unscaled),
    ('rm_width', float, 'synth', 'mom2CCFDF', u.rad/u.m**2),
    ('stokesI_model_coef', str, 'synth', 'polyCoeffs',None),

    # Less important quantities from the RMsynthesis (can be removed or modified after prototype verification?)
    ('noise_chan', float, 'synth', 'dQU', u.Jy/u.beam),
    ('fdf_noise_mad', float, 'synth', 'dFDFcorMAD', u.Jy/u.beam),
    ('fdf_noise_rms', float, 'synth', 'dFDFrms', u.Jy/u.beam),
    # We need to figure out what to do with these metrics, to simplify them.
    ('sigma_add_Q', float, 'synth', 'sigmaAddQ', u.dimensionless_unscaled),
    ('sigma_add_Q_err_plus', float, 'synth', 'dSigmaAddPlusQ', u.dimensionless_unscaled),
    ('sigma_add_Q_err_minus', float, 'synth','dSigmaAddMinusQ', u.dimensionless_unscaled),
    ('sigma_add_U', float, 'synth', 'sigmaAddU', u.dimensionless_unscaled),
    ('sigma_add_U_err_plus', float, 'synth','dSigmaAddPlusU', u.dimensionless_unscaled),
    ('sigma_add_U_err_minus', float, 'synth', 'dSigmaAddMinusU', u.dimensionless_unscaled),

    # Observation information. May be in header (what header?) or from meta data (what meta data, and how that gets pulled down, TBD)
    # these are position dependent in POSSUM: how to deal with?
    ('beam_maj', float, 'header', 'BMAJ', u.deg),
    ('beam_min', float, 'header', 'BMIN', u.deg),
    ('beam_pa', float, 'header', 'BPA', u.deg),

    # Not sure how these will be determined:
    # ('median_epoch_observation', float, '?', '?', cds.MJD),
    # ('integration_time', float, '?', '?', u.second),
    # ('interval_observation', float, '?', '?', u.day),

    # Metadata/flags?

    # Future:
    # Leakage estimate?
    # Stokes V estimates?
    # Faraday complexity flag
    # Data validation metrics
    # Metadata linking back to NRAO data (e.g, observation name)
]

# Jennifer's ideas for source validation metrics:
# Stokes I model fit failed
# unusually high noise
# number of channels flagged
# RM unusually high

output_cols = [x[0] for x in columns]
output_types = [x[1] for x in columns]
input_sources = [x[2] for x in columns]
input_names = [x[3] for x in columns]
output_units = [x[4] for x in columns]

# Columns in the expected input sourcelist (for RACS, BDSF)
sourcefinder_columns = (
                        # 'Gaussian_ID',
                        # 'Source_ID',
                        # 'Tile_ID',
                        # 'SBID',
                        # 'Obs_Start_Time',
                        # 'N_Gaus',
                        #'RA',
                        #'Dec',
                        #'E_RA',
                        #'E_Dec',
                        #'Total_flux_Gaussian',
                        #'E_Total_flux_Gaussian',
                        # 'Total_flux_Source',
                        # 'E_Total_flux_Source',
                        #'Peak_flux',
                        #'E_Peak_flux',
                        #'Maj',
                        #'E_Maj',
                        #'Min',
                        #'E_Min',
                        #'PA',
                        #'E_PA',
                        #'DC_Maj',
                        #'E_DC_Maj',
                        #'DC_Min',
                        #'E_DC_Min',
                        #'DC_PA',
                        #'E_DC_PA',
                        # 'S_Code',
                        # 'Separation_Tile_Centre',
                        # 'Noise',
                        #'Gal_lon',
                        #'Gal_lat'
                        )

# Values from Van Eck etal 
rmtab_column_descriptions = {
    "ra": "Right Ascension [ICRS]",
    "dec": "Declination [ICRS]",
    "l": "Galactic Longitude",
    "b": "Galactic Latitude",
    "pos_err": "Position uncertainty",
    "rm": "Rotation measure",
    "rm_err": "Error in RM",
    "rm_width": "Width in Faraday depth",
    "rm_width_err": "Error in width",
    "complex_flag": "Faraday complexity flag",
    "complex_test": "Faraday complexity metric",
    "rm_method": "RM determination method",
    "ionosphere": "Ionospheric correction method",
    "Ncomp": "Number of RM components",
    "stokesI": "Stokes I",
    "stokesI_err": "Error in Stokes I",
    "spectral_index": "Stokes I spectral index",
    "spectral_index_err": "Error in spectral index",
    "reffreq_I": "Reference frequency for Stokes I",
    "polint": "Polarized intensity",
    "polint_err": "Error in Pol.Int.",
    "pol_bias": "Polarization bias correction method",
    "flux_type": "Stokes extraction method",
    "aperture": "Integration aperture",
    "fracpol": "Fractional (linear) polarization",
    "fracpol_err": "Error in fractional polarization",
    "polangle": "Electric vector polarization angle",
    "polangle_err": "Error in EVPA",
    "reffreq_pol": "Reference frequency for polarization",
    "stokesQ": "Stokes Q",
    "stokesQ_err": "Error in Stokes Q",
    "stokesU": "Stokes U",
    "stokesU_err": "Error in Stokes U",
    "derot_polangle": "De-rotated EVPA",
    "derot_polangle_err": "Error in De-rotated EVPA",
    "stokesV": "Stokes V",
    "stokesV_err": "Error in Stokes V",
    "beam_maj": "Beam major axis",
    "beam_min": "Beam minor axis",
    "beam_pa": "Beam position angle",
    "reffreq_beam": "Reference frequency for beam",
    "minfreq": "Lowest frequency",
    "maxfreq": "Highest frequency",
    "channelwidth": "Typical channel width",
    "Nchan": "Number of channels",
    "rmsf_fwhm": "Full-width at half maximum of the RMSF",
    "noise_chan": "Typical per-channel noise in Q;U",
    "telescope": "Name of Telescope(s)",
    "int_time": "Integration time",
    "epoch": "Median epoch of observation",
    "interval": "Interval of observation",
    "leakage": "Instrumental leakage estimate",
    "beamdist": "Distance from beam centre - taken as centre of tile",
    "catalog": "Name of catalog",
    "dataref": "Data references",
    "cat_id": "Source ID in catalog",
    "type": "Source classification",
    "notes": "Notes",
}

extra_column_descriptions = {
    'ra_err': "Error in Right Ascension",
    'dec_err': "Error in Declination",
    'total_I_flux': "Total flux density in Stokes I",
    'total_I_flux_err': "Error in total flux density in Stokes I",
    'peak_I_flux': "Peak flux density in Stokes I",
    'peak_I_flux_err': "Error in peak flux density in Stokes I",
    'maj': "Major axis of Gaussian fit",
    'maj_err': "Error in major axis of Gaussian fit",
    'min': "Minor axis of Gaussian fit",
    'min_err': "Error in minor axis of Gaussian fit",
    'pa': "Position angle of Gaussian fit",
    'pa_err': "Error in position angle of Gaussian fit",
    'dc_maj': "Major axis of deconvolved Gaussian fit",
    'dc_maj_err': "Error in major axis of deconvolved Gaussian fit",
    'dc_min': "Minor axis of deconvolved Gaussian fit",
    'dc_min_err': "Error in minor axis of deconvolved Gaussian fit",
    'dc_pa': "Position angle of deconvolved Gaussian fit",
    'dc_pa_err': "Error in position angle of deconvolved Gaussian fit",
    'N_Gaus': "Number of Gaussians associated with source",
    'source_id': "Source ID",
    'tile_id': "Tile ID",
    'sbid': "Scheduling Block ID",
    'start_time': "Observation start time",
    'separation_tile_centre': "Separation from tile centre",
    's_code': "Source complexity classification",
    'fdf_noise_th': "Theoretical FDF noise",
    'refwave_sq_pol': "Reference wavelength squared",
    'snr_polint': "SNR of polarized intensity",
    'stokesI_model_coef': "Stokes I model coefficients",
    'fdf_noise_mad': "Median absolute deviation of FDF noise",
    'fdf_noise_rms': "RMS of FDF noise",
    'sigma_add': "Sigma_add complexity metric",
    'sigma_add_err': "Error in Sigma_add complexity metric",
    'snr_flag': "SNR flag",
    'leakage_flag': "Leakage flag",
    'channel_flag': "Channel flag",
    'stokes_I_fit_flag': "Stokes I fit flag",
    'complex_sigma_add_flag': "Sigma_add complexity flag",
    'complex_M2_CC_flag': "Second moment complexity flag",
}