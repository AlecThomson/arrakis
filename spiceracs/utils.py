#!/usr/bin/env python
"""Utility functions"""
import dataclasses
import functools
import json
import logging as log
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
import dask
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
from dask import delayed
from dask.delayed import Delayed
from distributed.client import futures_of
from distributed.diagnostics.progressbar import ProgressBar
from distributed.utils import LoopRunner, is_kernel
from FRion.correct import find_freq_axis
from scipy.optimize import curve_fit
from scipy.stats import normaltest
from spectral_cube import SpectralCube
from spectral_cube.utils import SpectralCubeWarning
from tornado.ioloop import IOLoop
from tqdm.auto import tqdm, trange

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)

print = functools.partial(print, flush=True)


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


def best_aic_func(aics, n_param):
    """Find the best AIC for a set of AICs using Occam's razor."""
    # Find the best AIC
    best_aic_idx = np.nanargmin(aics)
    best_aic = aics[best_aic_idx]
    best_n = n_param[best_aic_idx]
    log.debug(f"Lowest AIC is {best_aic}, with {best_n} params.")
    # Check if lower have diff < 2 in AIC
    aic_abs_diff = np.abs(aics - best_aic)
    bool_min_idx = np.zeros_like(aics).astype(bool)
    bool_min_idx[best_aic_idx] = True
    potential_idx = (aic_abs_diff[~bool_min_idx] < 2) & (
        n_param[~bool_min_idx] < best_n
    )
    if any(potential_idx):
        best_n = np.min(n_param[~bool_min_idx][potential_idx])
        best_aic_idx = np.where(n_param == best_n)[0][0]
        best_aic = aics[best_aic_idx]
        log.debug(
            f"Model within 2 of lowest AIC found. Occam says to take AIC of {best_aic}, with {best_n} params."
        )
    return best_aic, best_n, best_aic_idx


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
    x = nu / ref_nu
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
        aics = []
        params = []
        errors = []
        models = []
        highs = []
        lows = []
        fit_flags = []
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
                log.critical(f"Failed to fit {n}-term power law")
                aics.append(np.nan)
                params.append(np.ones_like(p0) * np.nan)
                errors.append(np.ones_like(p0) * np.nan)
                models.append(np.ones_like(freq))
                highs.append(np.ones_like(freq))
                lows.append(np.ones_like(freq))
                # 4 possible flags
                fit_flags.append(
                    {
                        "is_negative": True,
                        "is_not_finite": True,
                        "is_not_normal": True,
                        "is_close_to_zero": True,
                    }
                )
                continue
            best_p, covar = fit_res
            model_arr = model_func(freq, *best_p)
            model_high = model_func(freq, *(best_p + np.sqrt(np.diag(covar))))
            model_low = model_func(freq, *(best_p - np.sqrt(np.diag(covar))))
            model_err = model_high - model_low
            ssr = np.sum((flux[goodchan] - model_arr[goodchan]) ** 2)
            aic = akaike_info_criterion_lsq(ssr, len(p0), goodchan.sum())
            aics.append(aic)
            params.append(best_p)
            errors.append(np.sqrt(np.diag(covar)))
            models.append(model_arr)
            highs.append(model_high)
            lows.append(model_low)
            # Flag if model is negative
            is_negative = (model_arr < 0).any()
            if is_negative:
                log.warning(f"Stokes I flag: Model {n} is negative")
            # Flag if model is NaN or Inf
            is_not_finite = ~np.isfinite(model_arr).all()
            if is_not_finite:
                log.warning(f"Stokes I flag: Model {n} is not finite")
            # # Flag if model and data are statistically different
            residuals = flux[goodchan] - model_arr[goodchan]
            # Assume errors on resdiuals are the same as the data
            # i.e. assume the model has no error
            residuals_err = fluxerr[goodchan]
            residuals_norm = residuals / residuals_err
            # Test if the residuals are normally distributed
            ks,pval = normaltest(residuals_norm)
            is_not_normal = pval < 1e-6 # 1 in a million chance of being unlucky
            if is_not_normal:
                log.warning(f"Stokes I flag: Model {n} is not normally distributed - {pval=}, {ks=}")

            # Test if model is close to 0 within 1 sigma
            is_close_to_zero = (model_arr[goodchan] / fluxerr[goodchan] < 1).any()
            if is_close_to_zero:
                log.warning(f"Stokes I flag: Model {n} is close (1sigma) to 0")
            fit_flag = {
                    "is_negative": is_negative,
                    "is_not_finite": is_not_finite,
                    "is_not_normal": is_not_normal,
                    "is_close_to_zero": is_close_to_zero,
            }
            fit_flags.append(fit_flag)
            log.debug(f"{n}: {aic}")
        best_aic, best_n, best_aic_idx = best_aic_func(
            aics, np.array([n for n in range(nterms + 1)])
        )
        log.debug(f"Best fit: {best_n}, {best_aic}")
        best_p = params[best_n]
        best_e = errors[best_n]
        best_m = models[best_n]
        best_f = model_func_dict[best_n]
        best_flag = fit_flags[best_n]
        best_h = highs[best_n]
        best_l = lows[best_n]
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
        log.critical(f"Failed to fit power law: {e}")
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
                        "is_close_to_zero": True},
            ref_nu=np.nan,
            chi_sq=np.nan,
            chi_sq_red=np.nan,
        )


# stolen from https://stackoverflow.com/questions/32954486/zip-iterators-asserting-for-equal-length-in-python
def zip_equal(*iterables):
    sentinel = object()
    for combo in zip_longest(*iterables, fillvalue=sentinel):
        if sentinel in combo:
            raise ValueError("Iterables have different lengths")
        yield combo


def chunk_dask(
    outputs: list,
    client: distributed.Client,
    batch_size: int = 10_000,
    task_name="",
    progress_text="",
    verbose=True,
) -> list:
    chunk_outputs = []
    for i in trange(
        0, len(outputs), batch_size, desc=f"Chunking {task_name}", disable=(not verbose)
    ):
        outputs_chunk = outputs[i : i + batch_size]
        futures = client.persist(outputs_chunk)
        # dumb solution for https://github.com/dask/distributed/issues/4831
        if i == 0:
            log.debug("I sleep!")
            time.sleep(10)
            log.debug("I awake!")
        tqdm_dask(futures, desc=progress_text, disable=(not verbose))
        chunk_outputs.extend(futures)
    return chunk_outputs


def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """
    from math import sqrt

    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert columns in [1, 2]

    if fig_width is None:
        fig_width = 3.39 if columns == 1 else 6.9  # width in inches

    if fig_height is None:
        golden_mean = (sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
        fig_height = fig_width * golden_mean  # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print(
            "WARNING: fig_height too large:"
            + fig_height
            + "so will reduce to"
            + MAX_HEIGHT_INCHES
            + "inches."
        )
        fig_height = MAX_HEIGHT_INCHES

    params = {
        "backend": "pdf",
        "axes.labelsize": 8,  # fontsize for x and y labels (was 10)
        "axes.titlesize": 8,
        "font.size": 8,  # was 10
        "legend.fontsize": 8,  # was 10
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "text.usetex": False,
        "figure.figsize": [fig_width, fig_height],
        "font.family": "serif",
    }

    matplotlib.rcParams.update(params)


def delayed_to_da(list_of_delayed: List[Delayed], chunk: int = None) -> da.Array:
    """Convert list of delayed arrays to a dask array

    Args:
        list_of_delayed (List[delayed]): List of delayed objects
        chunk (int, optional): Chunksize to use. Defaults to None.

    Returns:
        da.Array: Dask array
    """
    sample = list_of_delayed[0].compute()
    dim = (len(list_of_delayed),) + sample.shape
    if chunk is None:
        c_dim = dim
    else:
        c_dim = (chunk,) + sample.shape
    darray_list = [
        da.from_delayed(lazy, dtype=sample.dtype, shape=sample.shape)
        for lazy in list_of_delayed
    ]
    darray = da.stack(darray_list, axis=0).reshape(dim).rechunk(c_dim)

    return darray


def yes_or_no(question: str) -> bool:
    """Ask a yes or no question via input()

    Args:
        question (str): Question to ask

    Returns:
        bool: True for yes, False for no
    """
    while "Please answer 'y' or 'n'":
        reply = str(input(question + " (y/n): ")).lower().strip()
        if reply[:1] == "y":
            ret = True
        if reply[:1] == "n":
            ret = False
    return ret


def fix_header(cutout_header: fits.Header, original_header: fits.Header) -> fits.Header:
    """Make cutout header the same as original header

    Args:
        cutout_header (fits.Header): Cutout header
        original_header (fits.Header): Original header

    Returns:
        fits.Header: Fixed header
    """
    axis_cut = find_freq_axis(cutout_header)
    axis_orig = find_freq_axis(original_header)
    fixed_header = cutout_header.copy()
    if axis_cut != axis_orig:
        for key, val in cutout_header.items():
            if key[-1] == str(axis_cut):
                fixed_header[f"{key[:-1]}{axis_orig}"] = val
                fixed_header[key] = original_header[key]

    return fixed_header


def deg_to_hms(deg: float) -> hms_tuple:
    """Convert degree to hms without astropy.

    Args:
        deg (float): Decimal degrees

    Returns:
        hms_tuple: HMS, like coord.ra.hms
    """
    h_per_d = 24 / 360
    hours = deg * h_per_d
    hour = float(int(hours))
    minutes = (hours - hour) * 60
    minute = float(int(minutes))
    seconds = (minutes - minute) * 60
    return hms_tuple(hour, minute, seconds)


def deg_to_dms(deg: float) -> dms_tuple:
    """Convert degree to hms without astropy.

    Args:
        deg (float): Decimal degrees

    Returns:
        hms_tuple: DMS, like coord.dec.dms
    """
    degree = float(int(deg))
    minutes = (deg - degree) * 60
    minute = float(int(minutes))
    seconds = (minutes - minute) * 60
    return dms_tuple(degree, minute, seconds)


def coord_to_string(coord: SkyCoord) -> Tuple[str, str]:
    """Convert coordinate to string without astropy

    Args:
        coord (SkyCoord): Coordinate

    Returns:
        Tuple[str,str]: Tuple of RA string, Dec string
    """
    ra = coord.ra
    dec = coord.dec

    ra_hms = deg_to_hms(ra.value)
    dec_dms = deg_to_dms(dec.value)

    ra_str = f"{ra_hms.h:02.0f}:{ra_hms.m:02.0f}:{ra_hms.s:06.3f}"
    dec_str = f"{dec_dms.d:02.0f}:{abs(dec_dms.m):02.0f}:{abs(dec_dms.s):05.2f}"
    return ra_str, dec_str


def test_db(
    host: str, username: str = None, password: str = None, verbose=True
) -> None:
    """Test connection to MongoDB

    Args:
        host (str): Mongo host IP.
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo password. Defaults to None.
        verbose (bool, optional): Verbose output. Defaults to True.

    Raises:
        Exception: If connection fails.
    """
    log.info("Testing MongoDB connection...")
    # default connection (ie, local)
    with pymongo.MongoClient(
        host=host,
        connect=False,
        username=username,
        password=password,
        authMechanism="SCRAM-SHA-256",
    ) as dbclient:  # type: pymongo.MongoClient
        try:
            dbclient.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")
        else:
            log.info("MongoDB connection succesful!")


def get_db(
    host: str, username: str = None, password: str = None
) -> Tuple[
    pymongo.collection.Collection,
    pymongo.collection.Collection,
    pymongo.collection.Collection,
]:
    """Get MongoDBs

    Args:
        host (str): Mongo host IP.
        username (str, optional): Username. Defaults to None.
        password (str, optional): Password. Defaults to None.

    Returns:
        Tuple[pymongo.Collection, pymongo.Collection, pymongo.Collection]: beams_col, island_col, comp_col
    """
    dbclient = pymongo.MongoClient(
        host=host,
        connect=False,
        username=username,
        password=password,
        authMechanism="SCRAM-SHA-256",
    )  # type: pymongo.MongoClient
    mydb = dbclient["spiceracs"]  # Create/open database
    comp_col = mydb["components"]  # Create/open collection
    island_col = mydb["islands"]  # Create/open collection
    beams_col = mydb["beams"]  # Create/open collection
    return beams_col, island_col, comp_col


def get_field_db(
    host: str, username=None, password=None
) -> pymongo.collection.Collection:
    """Get MongoDBs

    Args:
        host (str): Mongo host IP.
        username (str, optional): Username. Defaults to None.
        password (str, optional): Password. Defaults to None.

    Returns:
        pymongo.Collection: beams_col, island_col, comp_col
    """
    dbclient = pymongo.MongoClient(
        host=host,
        connect=False,
        username=username,
        password=password,
        authMechanism="SCRAM-SHA-256",
    )  # type: pymongo.MongoClient
    mydb = dbclient["spiceracs"]  # Create/open database
    field_col = mydb["fields"]  # Create/open collection
    return field_col


# stolen from https://github.com/tqdm/tqdm/issues/278
class TqdmProgressBar(ProgressBar):
    """Tqdm for Dask"""

    def __init__(
        self,
        keys,
        scheduler=None,
        interval="100ms",
        loop=None,
        complete=True,
        start=True,
        **tqdm_kwargs,
    ):
        super(TqdmProgressBar, self).__init__(keys, scheduler, interval, complete)
        self.tqdm = tqdm(keys, **tqdm_kwargs)
        self.loop = loop or IOLoop()

        if start:
            loop_runner = LoopRunner(self.loop)
            loop_runner.run_sync(self.listen)

    def _draw_bar(self, remaining, all, **kwargs):
        update_ct = (all - remaining) - self.tqdm.n
        self.tqdm.update(update_ct)

    def _draw_stop(self, **kwargs):
        self.tqdm.close()


def tqdm_dask(futures: distributed.Future, **kwargs) -> None:
    """Tqdm for Dask futures"""
    futures_list = futures_of(futures)
    if not isinstance(futures_list, (set, list)):
        futures_list = [futures_list]
    TqdmProgressBar(futures_list, **kwargs)


def port_forward(port: int, target: str) -> None:
    """Forward ports to local host

    Args:
        port (int): port to forward
        target (str): Target host
    """
    log.info(f"Forwarding {port} from localhost to {target}")
    cmd = f"ssh -N -f -R {port}:localhost:{port} {target}"
    command = shlex.split(cmd)
    output = subprocess.Popen(command)


def try_mkdir(dir_path: str, verbose=True):
    """Create directory if it doesn't exist

    Args:
        dir_path (str): Path to directory
        verbose (bool, optional): Verbose output. Defaults to True.
    """
    # Create output dir if it doesn't exist
    try:
        os.mkdir(dir_path)
        log.info(f"Made directory '{dir_path}'.")
    except FileExistsError:
        log.info(f"Directory '{dir_path}' exists.")


def try_symlink(src: str, dst: str, verbose=True):
    """Create symlink if it doesn't exist

    Args:
        src (str): Source path
        dst (str): Destination path
        verbose (bool, optional): Verbose output. Defaults to True.
    """
    # Create output dir if it doesn't exist
    try:
        os.symlink(src, dst)
        log.info(f"Made symlink '{dst}'.")
    except FileExistsError:
        log.info(f"Symlink '{dst}' exists.")


def head2dict(h: fits.Header) -> Dict[str, Any]:
    """Convert FITS header to a dict.

    Writes a cutout, as stored in source_dict, to disk. The file location
    should already be specified in source_dict. This format is intended
    for parallel use with pool.map syntax.

    Args:
        h: An astropy FITS header.

    Returns:
        data (dict): The FITS head converted to a dict.

    """
    data = {}
    for c in h.__dict__["_cards"]:
        if c[0] == "":
            continue
        data[c[0]] = c[1]
    return data


class MyEncoder(json.JSONEncoder):
    """Cutom JSON encorder.

    Parses the data stored in source_dict to JSON without
    errors.

    """

    def default(self, obj):  # pylint: disable=E0202
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.complex):
            return (obj.real, obj.imag)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, fits.Header):
            return head2dict(obj)
        elif dataclasses.is_dataclass(obj):
            return dataclasses.asdict(obj)
        else:
            return super(MyEncoder, self).default(obj)


def cpu_to_use(max_cpu: int, count: int) -> int:
    """Find number of cpus to use.

    Find the right number of cpus to use when dividing up a task, such
    that there are no remainders.

    Args:
        max_cpu (int): Maximum number of cores to use for a process.
        count (int): Number of tasks.

    Returns:
        Maximum number of cores to be used that divides into the number

    """
    factors = []
    for i in range(1, count + 1):
        if count % i == 0:
            factors.append(i)
    factors_arr = np.array(factors)
    return np.max(factors_arr[factors_arr <= max_cpu])


def getfreq(cube: str, outdir: str = None, filename: str = None):
    """Get list of frequencies from FITS data.

    Gets the frequency list from a given cube. Can optionally save
    frequency list to disk.

    Args:
        cube (str): File to get spectral axis from.

    Kwargs:
        outdir (str): Where to save the output file. If not given, data
            will not be saved to disk.

        filename (str): Name of frequency list file. Requires 'outdir'
            to also be specified.

        verbose (bool): Whether to print messages.

    Returns:
        freq (list): Frequencies of each channel in the input cube.

    """
    with fits.open(cube, memmap=True, mode="denywrite") as hdulist:
        hdu = hdulist[0]
        data = hdu.data
    wcs = WCS(hdu)
    freq = wcs.spectral.pixel_to_world(np.arange(data.shape[0]))  # Type: u.Quantity

    # Write to file if outdir is specified
    if outdir is None:
        return freq  # Type: u.Quantity
    else:
        if outdir[-1] == "/":
            outdir = outdir[:-1]
        if filename is None:
            outfile = f"{outdir}/frequencies.txt"
        else:
            outfile = f"{outdir}/{filename}"
        log.info(f"Saving to {outfile}")
        np.savetxt(outfile, np.array(freq))
        return freq, outfile  # Type: Tuple[u.Quantity, str]


def gettable(tabledir: str, keyword: str, verbose=True) -> Tuple[Table, str]:
    """Get a table from a directory given a keyword to glob.

    Args:
        tabledir (str): Directory.
        keyword (str): Keyword to glob for.
        verbose (bool, optional): Verbose output. Defaults to True.

    Returns:
        Tuple[Table, str]: Table and it's file location.
    """
    if tabledir[-1] == "/":
        tabledir = tabledir[:-1]
    # Glob out the necessary files
    files = glob(f"{tabledir}/*.{keyword}*.xml")  # Selvay VOTab
    filename = files[0]
    log.info(f"Getting table data from {filename}...")

    # Get selvay data from VOTab
    table = Table.read(filename, format="votable")
    table = table.to_pandas()
    str_df = table.select_dtypes([object])
    str_df = str_df.stack().str.decode("utf-8").unstack()
    for col in str_df:
        table[col] = str_df[col]
    return table, filename


def getdata(cubedir="./", tabledir="./", mapdata=None, verbose=True):
    """Get the spectral and source-finding data.

    Args:
        cubedir: Directory containing data cubes in FITS format.
        tabledir: Directory containing Selavy results.
        mapdata: 2D FITS image which corresponds to Selavy table.

    Kwargs:
        verbose (bool): Whether to print messages.

    Returns:
        datadict (dict): Dictionary of necessary astropy tables and
            Spectral cubes.

    """
    if cubedir[-1] == "/":
        cubedir = cubedir[:-1]

    if tabledir[-1] == "/":
        tabledir = tabledir[:-1]
    # Glob out the necessary files
    # Data cubes
    icubes = glob(f"{cubedir}/image.restored.i.*contcube*linmos.fits")
    qcubes = glob(f"{cubedir}/image.restored.q.*contcube*linmos.fits")
    ucubes = glob(f"{cubedir}/image.restored.u.*contcube*linmos.fits")
    vcubes = glob(f"{cubedir}/image.restored.v.*contcube*linmos.fits")

    cubes = [icubes, qcubes, ucubes, vcubes]
    # Selavy images
    selavyfits = mapdata
    # Get selvay data from VOTab
    i_tab, voisle = gettable(tabledir, "islands", verbose=verbose)  # Selvay VOTab
    components, tablename = gettable(tabledir, "components", verbose=verbose)

    log.info(f"Getting spectral data from: {cubes}\n")
    log.info(f"Getting source location data from: {selavyfits}\n")

    # Read data using Spectral cube
    i_taylor = SpectralCube.read(selavyfits, mode="denywrite")
    wcs_taylor = WCS(i_taylor.header)
    i_cube = SpectralCube.read(icubes[0], mode="denywrite")
    wcs_cube = WCS(i_cube.header)
    q_cube = SpectralCube.read(qcubes[0], mode="denywrite")
    u_cube = SpectralCube.read(ucubes[0], mode="denywrite")
    if len(vcubes) != 0:
        v_cube = SpectralCube.read(vcubes[0], mode="denywrite")
    else:
        v_cube = None
    # Mask out using Stokes I == 0 -- seems to be the current fill value
    mask = ~(i_cube == 0 * u.jansky / u.beam)
    i_cube = i_cube.with_mask(mask)
    mask = ~(q_cube == 0 * u.jansky / u.beam)
    q_cube = q_cube.with_mask(mask)
    mask = ~(u_cube == 0 * u.jansky / u.beam)
    u_cube = u_cube.with_mask(mask)

    datadict = {
        "i_tab": i_tab,
        "i_tab_comp": components,
        "i_taylor": i_taylor,
        "wcs_taylor": wcs_taylor,
        "wcs_cube": wcs_cube,
        "i_cube": i_cube,
        "q_cube": q_cube,
        "u_cube": u_cube,
        "v_cube": v_cube,
        "i_file": icubes[0],
        "q_file": qcubes[0],
        "u_file": ucubes[0],
        "v_file": vcubes[0],
    }

    return datadict


class Error(OSError):
    pass


class SameFileError(Error):
    """Raised when source and destination are the same file."""


class SpecialFileError(OSError):
    """Raised when trying to do a kind of operation (e.g. copying) which is
    not supported on a special file (e.g. a named pipe)"""


class ExecError(OSError):
    """Raised when a command could not be executed"""


class ReadError(OSError):
    """Raised when an archive cannot be read"""


class RegistryError(Exception):
    """Raised when a registry operation with the archiving
    and unpacking registeries fails"""


def _samefile(src, dst):
    # Macintosh, Unix.
    if hasattr(os.path, "samefile"):
        try:
            return os.path.samefile(src, dst)
        except OSError:
            return False


def copyfile(src, dst, *, follow_symlinks=True, verbose=True):
    """Copy data from src to dst.

    If follow_symlinks is not set and src is a symbolic link, a new
    symlink will be created instead of copying the file it points to.

    """
    if _samefile(src, dst):
        raise SameFileError("{!r} and {!r} are the same file".format(src, dst))

    for fn in [src, dst]:
        try:
            st = os.stat(fn)
        except OSError:
            # File most likely does not exist
            pass
        else:
            # XXX What about other special files? (sockets, devices...)
            if stat.S_ISFIFO(st.st_mode):
                raise SpecialFileError("`%s` is a named pipe" % fn)

    if not follow_symlinks and os.path.islink(src):
        os.symlink(os.readlink(src), dst)
    else:
        with open(src, "rb") as fsrc:
            with open(dst, "wb") as fdst:
                copyfileobj(fsrc, fdst, verbose=verbose)
    return dst


def copyfileobj(fsrc, fdst, length=16 * 1024, verbose=True):
    # copied = 0
    total = os.fstat(fsrc.fileno()).st_size
    with tqdm(
        total=total, disable=(not verbose), unit_scale=True, desc="Copying file"
    ) as pbar:
        while True:
            buf = fsrc.read(length)
            if not buf:
                break
            fdst.write(buf)
            copied = len(buf)
            pbar.update(copied)
