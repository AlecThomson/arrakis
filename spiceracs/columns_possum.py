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
    ('RA', float, 'cat', 'RA', u.deg),
    ('E_RA', float, 'cat', 'E_RA', u.arcsec),
    ('Dec', float, 'cat', 'Dec', u.deg),
    ('E_Dec', float, 'cat', 'E_Dec', u.arcsec),
    ('Total_I_flux_Gaussian', float, 'cat', 'Total_flux_Gaussian', u.mJy),
    ('E_Total_I_flux_Gaussian', float, 'cat',
     'E_Total_flux_Gaussian', u.mJy),
    ('Peak_I_flux_Gaussian', float, 'cat', 'Peak_flux', u.mJy/u.beam),
    ('E_Peak_I_flux_Gaussian', float, 'cat', 'E_Peak_flux', u.mJy/u.beam),
    # ('Gaussian_ID', str, 'cat', 'Gaussian_ID', None),
    ('Maj', float, 'cat', 'Maj', u.arcsec),
    ('E_Maj', float, 'cat', 'E_Maj', u.arcsec),
    ('Min', float, 'cat', 'Min', u.arcsec),
    ('E_Min', float, 'cat', 'E_Min', u.arcsec),
    ('PA', float, 'cat', 'PA', u.deg),
    ('E_PA', float, 'cat', 'E_PA', u.deg),
    ('DC_Maj',
     float, 'cat', 'DC_Maj', u.arcsec),
    ('E_DC_Maj',
     float, 'cat', 'E_DC_Maj', u.arcsec),
    ('DC_Min',
     float, 'cat', 'DC_Min', u.arcsec),
    ('E_DC_Min',
     float, 'cat', 'E_DC_Min', u.arcsec),
    ('DC_PA',
     float, 'cat', 'DC_PA', u.deg),
    ('E_DC_PA',
     float, 'cat', 'E_DC_PA', u.deg),
    ('Noise_I', float, 'cat', 'Noise', u.mJy/u.beam),
    ('Separation_Tile_Centre', float, 'cat', 'Separation_Tile_Centre', u.deg),

    # Important quantities from the RMsynthesis
    ('RM', float, 'synth', 'phiPeakPIfit_rm2', u.rad/u.m**2),
    ('E_RM', float, 'synth', 'dPhiPeakPIfit_rm2', u.rad/u.m**2),
    ('Polarised_intensity', float,
     'synth', 'ampPeakPIfitEff', u.Jy/u.beam),
    ('E_Polarised_intensity', float,
     'synth', 'dAmpPeakPIchan', u.Jy/u.beam),
    ('Stokes_Q', float, 'synth', 'peakFDFrealFit', u.Jy/u.beam),
    ('Stokes_U', float, 'synth', 'peakFDFimagFit', u.Jy/u.beam),
    # ('E_Stokes_QU', float, 'synth', 'dQU', u.Jy/u.beam),
    ('Polarisation_angle_0', float, 'synth', 'polAngle0Fit_deg', u.deg),
    ('E_Polarisation_angle_0', float, 'synth', 'dPolAngle0Fit_deg', u.deg),
    ('Fractional_polarisation', float, 'synth',
     'fracPol', u.dimensionless_unscaled),
    ('Reference_frequency', float, 'synth', 'freq0_Hz', u.Hz),
    ('Theoretical_FDF_noise', float, 'synth', 'dFDFth', u.Jy/u.beam),
    ('RMSF_FWHM', float, 'synth', 'fwhmRMSF', u.rad/u.m**2),
    ('Reference_wavelength_sq', float, 'synth', 'lam0Sq_m2', u.m**2),
    ('Reference_Stokes_I', float, 'synth', 'Ifreq0', u.Jy/u.beam),
    ('SNR_Polarised_intensity', float, 'synth', 'snrPIfit', u.dimensionless_unscaled),
    ('Lowest_frequency', float, 'synth', 'min_freq', u.Hz),
    ('Highest_frequency', float, 'synth', 'max_freq', u.Hz),
    ('Channel_width', float, 'synth', 'median_channel_width', u.Hz),
    ('Number_of_channels', int, 'synth',
     'N_channels', u.dimensionless_unscaled),
    ('M2_CC', float, 'synth', 'mom2CCFDF', u.rad/u.m**2),

    # Less important quantities from the RMsynthesis (can be removed or modified after prototype verification?)
    ('Typical_channel_noise', float, 'synth', 'dQU', u.Jy/u.beam),
    ('FDF_noise_MAD', float, 'synth', 'dFDFcorMAD', u.Jy/u.beam),
    ('FDF_noise_rms', float, 'synth', 'dFDFrms', u.Jy/u.beam),
    # We need to figure out what to do with these metrics, to simplify them.
    ('Sigma_add_Q', float, 'synth', 'sigmaAddQ', u.dimensionless_unscaled),
    ('E_Sigma_add_plus_Q', float, 'synth',
     'dSigmaAddPlusQ', u.dimensionless_unscaled),
    ('E_Sigma_add_minus_Q', float, 'synth',
     'dSigmaAddMinusQ', u.dimensionless_unscaled),
    ('Sigma_add_U', float, 'synth', 'sigmaAddU', u.dimensionless_unscaled),
    ('E_Sigma_add_plus_U', float, 'synth',
     'dSigmaAddPlusU', u.dimensionless_unscaled),
    ('E_Sigma_add_minus_U', float, 'synth',
     'dSigmaAddMinusU', u.dimensionless_unscaled),

    # Observation information. May be in header (what header?) or from meta data (what meta data, and how that gets pulled down, TBD)
    # these are position dependent in POSSUM: how to deal with?
    ('Beam_major_axis', float, 'header', 'BMAJ', u.deg),
    ('Beam_minor_axis', float, 'header', 'BMIN', u.deg),
    ('Beam_position_angle', float, 'header', 'BPA', u.deg),

    # Not sure how these will be determined:
    # ('median_epoch_observation', float, '?', '?', cds.MJD),
    # ('integration_time', float, '?', '?', u.second),
    # ('interval_observation', float, '?', '?', u.day),

    # Metadata/flags?
    ('Stokes_I_fit_flag', int, 'synth', 'IfitStat', u.dimensionless_unscaled),
    ('distance_from_beam_center', float, 'meta', 'beam_center_distance', u.deg)

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
sourcefinder_columns = ('Gaussian_ID',
                        'Source_ID',
                        'Tile_ID',
                        'SBID',
                        'Obs_Start_Time',
                        'N_Gaus',
                        #'RA',
                        #'Dec',
                        #'E_RA',
                        #'E_Dec',
                        #'Total_flux_Gaussian',
                        #'E_Total_flux_Gaussian',
                        'Total_flux_Source',
                        'E_Total_flux_Source',
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
                        'S_Code',
                        # 'Separation_Tile_Centre',
                        # 'Noise',
                        #'Gal_lon',
                        #'Gal_lat'
                        )
