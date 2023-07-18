#!/usr/bin/env python
"""Fitting utilities"""

import copy
import dataclasses
import functools
import json
import logging
import os
import shlex
import stat
import subprocess
import time
import warnings
from dataclasses import asdict, dataclass, make_dataclass
from functools import partial
from glob import glob
from itertools import zip_longest
from os import name
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import astropy.units as u
import dask.array as da
import dask.distributed as distributed
import numpy as np
import pymongo
from astropy.coordinates import SkyCoord
from astropy.coordinates.angles import dms_tuple, hms_tuple
from astropy.io import fits
from astropy.stats import akaike_info_criterion_lsq
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS
from casacore.tables import table
from casatasks import listobs
from dask import delayed
from dask.delayed import Delayed
from dask.distributed import Client, get_client
from distributed.client import futures_of
from distributed.diagnostics.progressbar import ProgressBar
from distributed.utils import LoopRunner, is_kernel
from FRion.correct import find_freq_axis
from prefect_dask import get_dask_client
from pymongo.collection import Collection
from scipy.optimize import curve_fit
from scipy.stats import normaltest
from spectral_cube import SpectralCube
from spectral_cube.utils import SpectralCubeWarning
from tornado.ioloop import IOLoop
from tqdm.auto import tqdm, trange

from arrakis.logger import logger

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)


def chi_squared(model: np.ndarray, data: np.ndarray, error: np.ndarray) -> float:
    """Calculate chi squared.

    Args:
        model (np.ndarray): Model flux.
        data (np.ndarray): Data flux.
        error (np.ndarray): Data error.

    Returns:
        np.ndarray: Chi squared.
    """
    return np.sum(((model - data) / error) ** 2)


def best_aic_func(aics: np.ndarray, n_param: np.ndarray) -> Tuple[float, int, int]:
    """Find the best AIC for a set of AICs using Occam's razor."""
    # Find the best AIC
    best_aic_idx = int(np.nanargmin(aics))
    best_aic = float(aics[best_aic_idx])
    best_n = int(n_param[best_aic_idx])
    logger.debug(f"Lowest AIC is {best_aic}, with {best_n} params.")
    # Check if lower have diff < 2 in AIC
    aic_abs_diff = np.abs(aics - best_aic)
    bool_min_idx = np.zeros_like(aics).astype(bool)
    bool_min_idx[best_aic_idx] = True
    potential_idx = (aic_abs_diff[~bool_min_idx] < 2) & (
        n_param[~bool_min_idx] < best_n
    )
    if not any(potential_idx):
        return best_aic, best_n, best_aic_idx

    bestest_n = int(np.min(n_param[~bool_min_idx][potential_idx]))
    bestest_aic_idx = int(np.where(n_param == bestest_n)[0][0])
    bestest_aic = float(aics[bestest_aic_idx])
    logger.debug(
        f"Model within 2 of lowest AIC found. Occam says to take AIC of {bestest_aic}, with {bestest_n} params."
    )
    return bestest_aic, bestest_n, bestest_aic_idx


# Stolen from GLEAM-X - thanks Uncle Timmy!
def power_law(nu: np.ndarray, norm: float, alpha: float, ref_nu: float) -> np.ndarray:
    """A power law model.

    Args:
        nu (np.ndarray): Frequency array.
        norm (float): Reference flux.
        alpha (float): Spectral index.
        ref_nu (float): Reference frequency.

    Returns:
        np.ndarray: Model flux.
    """
    return norm * (nu / ref_nu) ** alpha


def flat_power_law(nu: np.ndarray, norm: float, ref_nu: float) -> np.ndarray:
    """A flat power law model.

    Args:
        nu (np.ndarray): Frequency array.
        norm (float): Reference flux.
        ref_nu (float): Reference frequency.

    Returns:
        np.ndarray: Model flux.
    """
    x = ref_nu * np.ones_like(nu)
    return norm * x


def curved_power_law(
    nu: np.ndarray, norm: float, alpha: float, beta: float, ref_nu: float
) -> np.ndarray:
    """A curved power law model.

    Args:
        nu (np.ndarray): Frequency array.
        norm (float): Reference flux.
        alpha (float): Spectral index.
        beta (float): Spectral curvature.
        ref_nu (float): Reference frequency.

    Returns:
        np.ndarray: Model flux.
    """
    x = nu / ref_nu
    power = alpha + beta * np.log10(x)
    return norm * x**power


def fit_pl(
    freq: np.ndarray, flux: np.ndarray, fluxerr: np.ndarray, nterms: int
) -> dict:
    """Perform a power law fit to a spectrum.

    Args:
        freq (np.ndarray): Frequency array.
        flux (np.ndarray): Flux array.
        fluxerr (np.ndarray): Error array.
        nterms (int): Number of terms to use in the fit.

    Returns:
        dict: Best fit parameters.
    """
    try:
        goodchan = np.logical_and(
            np.isfinite(flux), np.isfinite(fluxerr)
        )  # Ignore NaN channels!
        ref_nu = np.nanmean(freq[goodchan])
        p0_long = (np.median(flux[goodchan]), -0.8, 0.0)
        model_func_dict = {
            0: partial(flat_power_law, ref_nu=ref_nu),
            1: partial(power_law, ref_nu=ref_nu),
            2: partial(curved_power_law, ref_nu=ref_nu),
        }

        # Initialise the save dict
        save_dict = {
            n: {} for n in range(nterms + 1)
        }  # type: Dict[int, Dict[str, Any]]
        for n in range(nterms + 1):
            p0 = p0_long[: n + 1]
            save_dict[n]["aics"] = np.nan
            save_dict[n]["params"] = np.ones_like(p0) * np.nan
            save_dict[n]["errors"] = np.ones_like(p0) * np.nan
            save_dict[n]["models"] = np.ones_like(freq)
            save_dict[n]["highs"] = np.ones_like(freq)
            save_dict[n]["lows"] = np.ones_like(freq)
            # 4 possible flags
            save_dict[n]["fit_flags"] = {
                "is_negative": True,
                "is_not_finite": True,
                "is_not_normal": True,
                "is_close_to_zero": True,
            }

        # Now iterate over the number of terms
        for n in range(nterms + 1):
            p0 = p0_long[: n + 1]
            model_func = model_func_dict[n]
            try:
                fit_res = curve_fit(
                    model_func,
                    freq[goodchan],
                    flux[goodchan],
                    p0=p0,
                    sigma=fluxerr[goodchan],
                    absolute_sigma=True,
                )
            except RuntimeError:
                logger.critical(f"Failed to fit {n}-term power law")
                continue

            best, covar = fit_res
            model_arr = model_func(freq, *best)
            model_high = model_func(freq, *(best + np.sqrt(np.diag(covar))))
            model_low = model_func(freq, *(best - np.sqrt(np.diag(covar))))
            model_err = model_high - model_low
            ssr = np.sum((flux[goodchan] - model_arr[goodchan]) ** 2)
            aic = akaike_info_criterion_lsq(ssr, len(p0), goodchan.sum())

            # Save the results
            save_dict[n]["aics"] = aic
            save_dict[n]["params"] = best
            save_dict[n]["errors"] = np.sqrt(np.diag(covar))
            save_dict[n]["models"] = model_arr
            save_dict[n]["highs"] = model_high
            save_dict[n]["lows"] = model_low

            # Calculate the flags
            # Flag if model is negative
            is_negative = (model_arr < 0).any()
            if is_negative:
                logger.warning(f"Stokes I flag: Model {n} is negative")
            # Flag if model is NaN or Inf
            is_not_finite = ~np.isfinite(model_arr).all()
            if is_not_finite:
                logger.warning(f"Stokes I flag: Model {n} is not finite")
            # # Flag if model and data are statistically different
            residuals = flux[goodchan] - model_arr[goodchan]
            # Assume errors on resdiuals are the same as the data
            # i.e. assume the model has no error
            residuals_err = fluxerr[goodchan]
            residuals_norm = residuals / residuals_err
            # Test if the residuals are normally distributed
            ks, pval = normaltest(residuals_norm)
            is_not_normal = pval < 1e-6  # 1 in a million chance of being unlucky
            if is_not_normal:
                logger.warning(
                    f"Stokes I flag: Model {n} is not normally distributed - {pval=}, {ks=}"
                )

            # Test if model is close to 0 within 1 sigma
            is_close_to_zero = (model_arr[goodchan] / fluxerr[goodchan] < 1).any()
            if is_close_to_zero:
                logger.warning(f"Stokes I flag: Model {n} is close (1sigma) to 0")
            fit_flag = {
                "is_negative": is_negative,
                "is_not_finite": is_not_finite,
                "is_not_normal": is_not_normal,
                "is_close_to_zero": is_close_to_zero,
            }
            save_dict[n]["fit_flags"] = fit_flag
            logger.debug(f"{n}: {aic}")

        # Now find the best model
        best_aic, best_n, best_aic_idx = best_aic_func(
            np.array([save_dict[n]["aics"] for n in range(nterms + 1)]),
            np.array([n for n in range(nterms + 1)]),
        )
        logger.debug(f"Best fit: {best_n}, {best_aic}")
        best_p = save_dict[best_n]["params"]
        best_e = save_dict[best_n]["errors"]
        best_m = save_dict[best_n]["models"]
        best_f = model_func_dict[best_n]
        best_flag = save_dict[best_n]["fit_flags"]
        best_h = save_dict[best_n]["highs"]
        best_l = save_dict[best_n]["lows"]
        chi_sq = chi_squared(
            model=best_m[goodchan],
            data=flux[goodchan],
            error=fluxerr[goodchan],
        )
        chi_sq_red = chi_sq / (goodchan.sum() - len(best_p))
        return dict(
            best_n=best_n,
            best_p=best_p,
            best_e=best_e,
            best_m=best_m,
            best_h=best_h,
            best_l=best_l,
            best_f=best_f,
            fit_flag=best_flag,
            ref_nu=ref_nu,
            chi_sq=chi_sq,
            chi_sq_red=chi_sq_red,
        )
    except Exception as e:
        logger.critical(f"Failed to fit power law: {e}")
        return dict(
            best_n=np.nan,
            best_p=[np.nan],
            best_e=[np.nan],
            best_m=np.ones_like(freq),
            best_h=np.ones_like(freq),
            best_l=np.ones_like(freq),
            best_f=None,
            fit_flag={
                "is_negative": True,
                "is_not_finite": True,
                "is_not_normal": True,
                "is_close_to_zero": True,
            },
            ref_nu=np.nan,
            chi_sq=np.nan,
            chi_sq_red=np.nan,
        )
