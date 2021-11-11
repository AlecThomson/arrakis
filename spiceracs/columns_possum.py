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
    ('e_ra', float, 'cat', 'E_RA', u.arcsec),
    ('dec', float, 'cat', 'Dec', u.deg),
    ('e_dec', float, 'cat', 'E_Dec', u.arcsec),
    ('total_flux_sourcefinder', float, 'cat', 'Total_flux_Gaussian', 1e-3*u.Jy),
    ('e_total_flux_sourcefinder', float, 'cat',
     'E_Total_flux_Gaussian', 1e-3*u.Jy),
    ('peak_flux_sourcefinder', float, 'cat', 'Peak_flux', 1e-3*u.Jy/u.beam),
    ('e_peak_flux_sourcefinder', float, 'cat', 'E_Peak_flux', 1e-3*u.Jy/u.beam),
    ('id_number_sourcefinder', str, 'cat', 'Gaussian_ID', None),
    ('major_axis_sourcefinder', float, 'cat', 'Maj', u.arcsec),
    ('e_major_axis_sourcefinder', float, 'cat', 'E_Maj', u.arcsec),
    ('minor_axis_sourcefinder', float, 'cat', 'Min', u.arcsec),
    ('e_minor_axis_sourcefinder', float, 'cat', 'E_Min', u.arcsec),
    ('position_angle_sourcefinder', float, 'cat', 'PA', u.deg),
    ('e_position_angle_sourcefinder', float, 'cat', 'E_PA', u.deg),
    ('deconvolved_major_axis_sourcefinder',
     float, 'cat', 'DC_Maj', u.arcsec),
    ('e_deconvolved_major_axis_sourcefinder',
     float, 'cat', 'E_DC_Maj', u.arcsec),
    ('deconvolved_minor_axis_sourcefinder',
     float, 'cat', 'DC_Min', u.arcsec),
    ('e_deconvolved_minor_axis_sourcefinder',
     float, 'cat', 'E_DC_Min', u.arcsec),
    ('deconvolved_position_angle_sourcefinder',
     float, 'cat', 'DC_PA', u.deg),
    ('e_deconvolved_position_angle_sourcefinder',
     float, 'cat', 'E_DC_PA', u.deg),
    ('local_I_rms_sourcefinder', float, 'cat', 'Noise', 1e-3*u.Jy/u.beam),

    # Important quantities from the RMsynthesis
    ('peak_RM', float, 'synth', 'phiPeakPIfit_rm2', u.rad/u.m**2),
    ('e_peak_RM', float, 'synth', 'dPhiPeakPIfit_rm2', u.rad/u.m**2),
    ('peak_polarized_intensity', float,
     'synth', 'ampPeakPIfitEff', u.Jy/u.beam),
    ('e_peak_polarized_intensity', float,
     'synth', 'dAmpPeakPIchan', u.Jy/u.beam),
    ('peak_Q', float, 'synth', 'peakFDFrealFit', u.Jy/u.beam),
    ('peak_U', float, 'synth', 'peakFDFimagFit', u.Jy/u.beam),
    ('err_QU', float, 'synth', 'dQU', u.Jy/u.beam),
    ('derotated_angle', float, 'synth', 'polAngle0Fit_deg', u.deg),
    ('e_derotated_angle', float, 'synth', 'dPolAngle0Fit_deg', u.deg),
    ('fractional_polarization', float, 'synth',
     'fracPol', u.dimensionless_unscaled),
    ('reference_frequency_I_fracpol', float, 'synth', 'freq0_Hz', u.Hz),
    ('theoretical_noise', float, 'synth', 'dFDFth', u.Jy/u.beam),
    ('rmsf_fwhm', float, 'synth', 'fwhmRMSF', u.rad/u.m**2),
    ('lam0Sq_m2', float, 'synth', 'lam0Sq_m2', u.m**2),
    ('total_flux_fit', float, 'synth', 'Ifreq0', u.Jy/u.beam),
    ('snr_fit', float, 'synth', 'snrPIfit', u.dimensionless_unscaled),
    ('lowest_frequency', float, 'synth', 'min_freq', u.Hz),
    ('highest_frequency', float, 'synth', 'max_freq', u.Hz),
    ('channel_width', float, 'synth', 'median_channel_width', u.Hz),
    ('number_of_channels', int, 'synth',
     'N_channels', u.dimensionless_unscaled),

    # Less important quantities from the RMsynthesis (can be removed or modified after prototype verification?)
    ('typical_channel_noise', float, 'synth', 'dQU', u.Jy/u.beam),
    ('spectrum_noise_corMAD', float, 'synth', 'dFDFcorMAD', u.Jy/u.beam),
    ('spectrum_noise_rms', float, 'synth', 'dFDFrms', u.Jy/u.beam),
    # We need to figure out what to do with these metrics, to simplify them.
    ('sigmaAddQ', float, 'synth', 'sigmaAddQ', u.dimensionless_unscaled),
    ('dSigmaAddPlusQ', float, 'synth',
     'dSigmaAddPlusQ', u.dimensionless_unscaled),
    ('dSigmaAddMinusQ', float, 'synth',
     'dSigmaAddMinusQ', u.dimensionless_unscaled),
    ('sigmaAddU', float, 'synth', 'sigmaAddU', u.dimensionless_unscaled),
    ('dSigmaAddPlusU', float, 'synth',
     'dSigmaAddPlusU', u.dimensionless_unscaled),
    ('dSigmaAddMinusU', float, 'synth',
     'dSigmaAddMinusU', u.dimensionless_unscaled),

    # Observation information. May be in header (what header?) or from meta data (what meta data, and how that gets pulled down, TBD)
    # these are position dependent in POSSUM: how to deal with?
    ('beam_major_axis', float, 'header', 'BMAJ', u.deg),
    ('beam_minor_axis', float, 'header', 'BMIN', u.deg),
    ('beam_position_angle', float, 'header', 'BPA', u.deg),

    # Not sure how these will be determined:
    ('median_epoch_observation', float, '?', '?', cds.MJD),
    ('integration_time', float, '?', '?', u.second),
    ('interval_observation', float, '?', '?', u.day),

    # Metadata/flags?
    ('stokesI_fit_flag', int, 'synth', 'IfitStat', u.dimensionless_unscaled),
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
                        'Separation_Tile_Centre',
                        'Noise',
                        #'Gal_lon',
                        #'Gal_lat'
                        )
