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
    ('reffreq_I', float, 'synth', 'poly_reffreq', u.Hz),
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
    ('stokesI_model_coef_err', str, 'synth', 'polyCoefferr', None),
    ('stokesI_model_order', int, 'synth', 'polyOrd', u.dimensionless_unscaled),

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

extra_column_descriptions = {
    'ra_err': {
        "description":"Error in Right Ascension",
        "ucd": "stat.error;pos.eq.ra",
    },
    'dec_err': {
        "description":"Error in Declination",
        "ucd": "stat.error;pos.eq.dec",
    },
    'total_I_flux': {
        "description":"Total flux density in Stokes I",
        "ucd": "phot.flux.density;arith.sum;phys.polarization.stokes.I",
    },
    'total_I_flux_err': {
        "description":"Error in total flux density in Stokes I",
        "ucd": "stat.error;phot.flux.density;arith.sum;phys.polarization.stokes.I",
    },
    'peak_I_flux': {
        "description":"Peak flux density in Stokes I",
        "ucd": "phot.flux.density;stat.max;phys.polarization.stokes.I",
    },
    'peak_I_flux_err': {
        "description":"Error in peak flux density in Stokes I",
        "ucd": "stat.error;phot.flux.density;stat.max;phys.polarization.stokes.I",
    },
    'maj': {
        "description":"Major axis of Gaussian fit",
        "ucd": "phys.angSize.smajAxis",
    },
    'maj_err': {
        "description":"Error in major axis of Gaussian fit",
        "ucd": "stat.error;phys.angSize.smajAxis",
    },
    'min': {
        "description":"Minor axis of Gaussian fit",
        "ucd": "phys.angSize.sminAxis",
    },
    'min_err': {
        "description":"Error in minor axis of Gaussian fit",
        "ucd": "stat.error;phys.angSize.sminAxis",
    },
    'pa': {
        "description":"Position angle of Gaussian fit",
        "ucd": "phys.angSize;pos.posAng",
    },
    'pa_err': {
        "description":"Error in position angle of Gaussian fit",
        "ucd": "stat.error;phys.angSize;pos.posAng",
    },
    'dc_maj': {
        "description":"Major axis of deconvolved Gaussian fit",
        "ucd": "phys.angSize.smajAxis",
    },
    'dc_maj_err': {
        "description":"Error in major axis of deconvolved Gaussian fit",
        "ucd": "stat.error;phys.angSize.smajAxis",
    },
    'dc_min': {
        "description":"Minor axis of deconvolved Gaussian fit",
        "ucd": "phys.angSize.sminAxis",
    },
    'dc_min_err': {
        "description":"Error in minor axis of deconvolved Gaussian fit",
        "ucd": "stat.error;phys.angSize.sminAxis",
    },
    'dc_pa': {
        "description":"Position angle of deconvolved Gaussian fit",
        "ucd": "phys.angSize;pos.posAng",
    },
    'dc_pa_err': {
        "description":"Error in position angle of deconvolved Gaussian fit",
        "ucd": "stat.error;phys.angSize;pos.posAng",
    },
    'N_Gaus': {
        "description":"Number of Gaussians associated with source",
        "ucd": "meta.number",
    },
    'source_id': {
        "description":"Source ID",
        "ucd": "meta.id",
    },
    'tile_id': {
        "description":"Tile ID",
        "ucd": "meta.id",
    },
    'sbid': {
        "description":"Scheduling Block ID",
        "ucd": "meta.id",
    },
    'start_time': {
        "description":"Observation start time",
        "ucd": "time.start",
    },
    'separation_tile_centre': {
        "description":"Separation from tile centre",
        "ucd": "pos.angDistance",
    },
    's_code': {
        "description":"Source complexity classification",
        "ucd": "meta.code.class",
    },
    'fdf_noise_th': {
        "description":"Theoretical FDF noise",
        "ucd": "stat.error;phot.flux.density;phys.polarization.linear",
    },
    'refwave_sq_pol': {
        "description":"Reference wavelength squared",
        "ucd": "em.wl;arith.squared",
    },
    'snr_polint': {
        "description":"SNR of polarized intensity",
        "ucd": "stat.snr;phot.flux.density;phys.polarization.linear",
    },
    'stokesI_model_coef': {
        "description":"Stokes I model coefficients",
        "ucd": "stat.fit.param;phys.polarization.stokes.I",
    },
    'fdf_noise_mad': {
        "description":"Median absolute deviation of FDF noise",
        "ucd": "stat.error;phot.flux.density;phys.polarization.linear",
    },
    'fdf_noise_rms': {
        "description":"RMS of FDF noise",
        "ucd": "stat.error;phot.flux.density;phys.polarization.linear",
    },
    'sigma_add': {
        "description":"Sigma_add complexity metric",
        "ucd": "phys.polarization.rotMeasure;stat.stdev",
    },
    'sigma_add_err': {
        "description":"Error in Sigma_add complexity metric",
        "ucd": "stat.error;phys.polarization.rotMeasure;stat.stdev",
    },
    'snr_flag': {
        "description":"SNR flag",
        "ucd": "meta.code.qual",
    },
    'leakage_flag': {
        "description":"Leakage flag",
        "ucd": "meta.code.qual",
    },
    'channel_flag': {
        "description":"Channel flag",
        "ucd": "meta.code.qual",
    },
    'stokes_I_fit_flag': {
        "description":"Stokes I fit flag",
        "ucd": "meta.code.qual",
    },
    'complex_sigma_add_flag': {
        "description":"Sigma_add complexity flag",
        "ucd": "meta.code",
    },
    'complex_M2_CC_flag': {
        "description":"Second moment complexity flag",
        "ucd": "meta.code",
    },
}