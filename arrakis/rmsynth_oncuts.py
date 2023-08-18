#!/usr/bin/env python3
"""Run RM-CLEAN on cutouts in parallel"""
import logging
import os
import traceback
import warnings
from glob import glob
from pathlib import Path
from pprint import pformat
from shutil import copyfile
from typing import List
from typing import NamedTuple as Struct
from typing import Optional, Tuple, Union

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymongo
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling import models
from astropy.stats import mad_std, sigma_clip
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from dask import delayed
from dask.distributed import Client, LocalCluster
from IPython import embed
from radio_beam import Beam
from RMtools_1D import do_RMsynth_1D
from RMtools_3D import do_RMsynth_3D
from RMutils.util_misc import create_frac_spectra
from scipy.stats import norm
from tqdm import tqdm, trange

from arrakis.logger import logger
from arrakis.utils.database import get_db, test_db
from arrakis.utils.fitsutils import getfreq
from arrakis.utils.fitting import fit_pl
from arrakis.utils.io import try_mkdir
from arrakis.utils.pipeline import chunk_dask, logo_str

logger.setLevel(logging.INFO)


class Spectrum(Struct):
    """Single spectrum"""

    data: np.ndarray
    """The spectrum data"""
    rms: np.ndarray
    """The RMS of the spectrum"""
    bkg: np.ndarray
    """The background of the spectrum"""
    filename: str
    """The filename associated with the spectrum"""
    header: fits.Header
    """The header associated with the spectrum"""


class StokesSpectra(Struct):
    """Multi Stokes spectra"""

    i: Spectrum
    """The Stokes I spectrum"""
    q: Spectrum
    """The Stokes Q spectrum"""
    u: Spectrum
    """The Stokes U spectrum"""


class StokesIFitResult(Struct):
    """Stokes I fit results"""

    alpha: Optional[float]
    """The alpha parameter of the fit"""
    amplitude: Optional[float]
    """The amplitude parameter of the fit"""
    x_0: Optional[float]
    """The x_0 parameter of the fit"""
    model_repr: Optional[str]
    """The model representation of the fit"""
    modStokesI: Optional[np.ndarray]
    """The model Stokes I spectrum"""
    fit_dict: Optional[dict]
    """The dictionary of the fit results"""


@delayed
def rmsynthoncut3d(
    island_id: str,
    beam: dict,
    outdir: str,
    freq: np.ndarray,
    field: str,
    phiMax_radm2: Optional[float] = None,
    dPhi_radm2: Optional[float] = None,
    nSamples: int = 5,
    weightType: str = "variance",
    fitRMSF: bool = True,
    not_RMSF: bool = False,
    rm_verbose: bool = False,
    ion: bool = False,
):
    """3D RM-synthesis

    Args:
        island_id (str): RACS Island ID
        freq (list): Frequencies in Hz
        host (str): Host of MongoDB
        field (str): RACS field ID
        database (bool, optional): Update MongoDB. Defaults to False.
        phiMax_radm2 (float, optional): Max Faraday depth. Defaults to None.
        dPhi_radm2 (float, optional): Faraday dpeth channel width. Defaults to None.
        nSamples (int, optional): Samples acorss RMSF. Defaults to 5.
        weightType (str, optional): Weighting type. Defaults to 'variance'.
        fitRMSF (bool, optional): Fit RMSF. Defaults to False.
        not_RMSF (bool, optional): Skip calculation of RMSF. Defaults to False.
        rm_verbose (bool, optional): Verbose RMsynth. Defaults to False.
    """

    iname = island_id
    ifile = os.path.join(outdir, beam["beams"][field]["i_file"])

    if ion:
        qfile = os.path.join(outdir, beam["beams"][field]["q_file_ion"])
        ufile = os.path.join(outdir, beam["beams"][field]["u_file_ion"])
    else:
        qfile = os.path.join(outdir, beam["beams"][field]["q_file"])
        ufile = os.path.join(outdir, beam["beams"][field]["u_file"])
    # vfile = beam['beams'][field]['v_file']

    header, dataQ = do_RMsynth_3D.readFitsCube(qfile, rm_verbose)
    header, dataU = do_RMsynth_3D.readFitsCube(ufile, rm_verbose)
    header, dataI = do_RMsynth_3D.readFitsCube(ifile, rm_verbose)

    dataQ = np.squeeze(dataQ)
    dataU = np.squeeze(dataU)
    dataI = np.squeeze(dataI)

    if np.isnan(dataI).all() or np.isnan(dataQ).all() or np.isnan(dataU).all():
        logger.critical(f"Cubelet {iname} is entirely NaN")
        myquery = {"Source_ID": iname}
        badvalues = {
            "$set": {
                "rmsynth3d": False,
            }
        }
        return pymongo.UpdateOne(myquery, badvalues)

    bkgq, rmsq = cubelet_bane(dataQ, header)
    rmsq[rmsq == 0] = np.nan
    rmsq[np.isnan(rmsq)] = np.nanmedian(rmsq)

    bkgu, rmsu = cubelet_bane(dataU, header)
    rmsu[rmsu == 0] = np.nan
    rmsu[np.isnan(rmsu)] = np.nanmedian(rmsu)
    rmsArr = np.max([rmsq, rmsu], axis=0)

    # Run 3D RM-synthesis on the cubes
    dataArr = do_RMsynth_3D.run_rmsynth(
        dataQ=dataQ,
        dataU=dataU,
        freqArr_Hz=freq,
        dataI=None,
        rmsArr=rmsArr,
        phiMax_radm2=phiMax_radm2,
        dPhi_radm2=dPhi_radm2,
        nSamples=nSamples,
        weightType=weightType,
        fitRMSF=fitRMSF,
        nBits=32,
        verbose=rm_verbose,
        not_rmsf=not_RMSF,
    )

    prefix = f"{iname}_"
    # Write to files
    do_RMsynth_3D.writefits(
        dataArr,
        headtemplate=header,
        fitRMSF=fitRMSF,
        prefixOut=prefix,
        outDir=os.path.dirname(ifile),
        write_seperate_FDF=True,
        not_rmsf=not_RMSF,
        nBits=32,
        verbose=rm_verbose,
    )

    myquery = {"Source_ID": iname}
    # Prep header
    head_dict = dict(header)
    head_dict.pop("", None)
    if "COMMENT" in head_dict.keys():
        head_dict["COMMENT"] = str(head_dict["COMMENT"])

    outer_dir = os.path.basename(os.path.dirname(ifile))

    newvalues = {
        "$set": {
            "rm3dfiles": {
                "FDF_real_dirty": os.path.join(
                    outer_dir, f"{prefix}FDF_real_dirty.fits"
                ),
                "FDF_im_dirty": os.path.join(outer_dir, f"{prefix}FDF_im_dirty.fits"),
                "FDF_tot_dirty": os.path.join(outer_dir, f"{prefix}FDF_tot_dirty.fits"),
                "RMSF_real": os.path.join(outer_dir, f"{prefix}RMSF_real.fits"),
                "RMSF_tot": os.path.join(outer_dir, f"{prefix}RMSF_tot.fits"),
                "RMSF_FWHM": os.path.join(outer_dir, f"{prefix}RMSF_FWHM.fits"),
            },
            "rmsynth3d": True,
            "header": dict(header),
        }
    }
    return pymongo.UpdateOne(myquery, newvalues)


def cubelet_bane(cubelet: np.ndarray, header: fits.Header) -> Tuple[np.ndarray]:
    """Background and noise estimation on a cubelet

    Args:
        cubelet (np.ndarray): 3D array of data
        header (fits.Header): Header of cubelet

    Returns:
        Tuple[np.ndarray]: Background and noise per channel
    """
    # Get beam and pixel information
    wcs = WCS(header).celestial
    pixelscales = max(proj_plane_pixel_scales(wcs)) * u.deg / u.pixel
    beam = Beam.from_fits_header(header)
    pix_per_beam = (beam.minor / pixelscales).decompose()

    # Create annulus mask
    x_grid, y_grid = np.meshgrid(
        np.arange(cubelet.shape[2]), np.arange(cubelet.shape[1])
    )
    x_grid -= cubelet.shape[2] // 2
    y_grid -= cubelet.shape[1] // 2
    r_grid = np.hypot(x_grid, y_grid)
    mask = (r_grid > 1.5 * pix_per_beam.to(u.pix).value) & (
        r_grid < 5 * pix_per_beam.to(u.pix).value
    )
    data_masked = cubelet[:, mask]

    # Fit background and noise for each channel
    background = np.zeros(cubelet.shape[0]) * np.nan
    noise = np.zeros(cubelet.shape[0]) * np.nan
    for chan, plane in enumerate(data_masked):
        plane = plane[np.isfinite(plane)]
        if len(plane) == 0:
            continue
        background[chan], noise[chan] = norm.fit(plane)

    return background, noise


# TODO: Add NxN sum of pixels
def extract_single_spectrum(
    coord: SkyCoord,
    stokes: str,
    ion: bool,
    field_dict: dict,
    outdir: str,
) -> Spectrum:
    """Extract a single spectrum from a cubelet"""
    if ion and (stokes == "q" or stokes == "u"):
        key = f"{stokes}_file_ion"
    else:
        key = f"{stokes}_file"
    filename = os.path.join(outdir, field_dict[key])
    with fits.open(filename, mode="denywrite", memmap=True) as hdulist:
        hdu = hdulist[0]
        data = np.squeeze(hdu.data)
        header = hdu.header

    bkg, rms = cubelet_bane(data, header)
    rms[np.isnan(rms)] = np.nanmedian(rms)
    wcs = WCS(header)
    x, y = np.array(wcs.celestial.world_to_pixel(coord)).round().astype(int)
    spectrum_arr = data[:, y, x]
    spectrum_arr[spectrum_arr == 0] = np.nan
    # Do background subtraction
    spectrum_arr -= bkg
    return Spectrum(
        data=spectrum_arr,
        rms=rms,
        bkg=bkg,
        filename=filename,
        header=header,
    )


def extract_all_spectra(
    coord: SkyCoord,
    ion: bool,
    field_dict: dict,
    outdir: str,
) -> StokesSpectra:
    """Extract spectra from cubelets"""
    StokesSpectra(
        *[
            extract_single_spectrum(
                coord=coord,
                stokes=stokes,
                ion=ion,
                field_dict=field_dict,
                outdir=outdir,
            )
            for stokes in "iqu"
        ]
    )


def sigma_clip_spectra(
    stokes_spectra: StokesSpectra,
) -> StokesSpectra:
    """Sigma clip spectra

    Find outliers in the RMS spectra and set them to NaN

    Args:
        stokes_spectra (StokesSpectra): The Stokes spectra

    Returns:
        StokesSpectra: The filtered Stokes spectra
    """
    filter_list: List[np.ndarry] = []
    for spectrum in stokes_spectra:
        rms_filter = sigma_clip(
            spectrum.rms,
            sigma=5,
            stdfunc=mad_std,
            cenfunc=np.nanmedian,
        )
        filter_list.append(rms_filter.mask)
    filter_idx = np.any(filter_list, axis=0)

    filtered_data_list: List[Spectrum] = []
    for spectrum in stokes_spectra:
        filtered_data = spectrum.data.copy()
        filtered_data[filter_idx] = np.nan
        filtered_spectrum = Spectrum(
            data=filtered_data,
            rms=spectrum.rms,
            bkg=spectrum.bkg,
            filename=spectrum.filename,
            header=spectrum.header,
        )
        filtered_data_list.append(filtered_spectrum)
    return StokesSpectra(*filtered_data_list)


def fit_stokes_I(
    freq: np.ndarray,
    coord: SkyCoord,
    tt0: Optional[str] = None,
    tt1: Optional[str] = None,
    do_own_fit: bool = False,
    iarr: Optional[np.ndarray] = None,
    rmsi: Optional[np.ndarray] = None,
    polyOrd: Optional[int] = None,
) -> StokesIFitResult:
    if tt0 and tt1:
        mfs_i_0 = fits.getdata(tt0, memmap=True)
        mfs_i_1 = fits.getdata(tt1, memmap=True)
        mfs_head = fits.getheader(tt0)
        mfs_wcs = WCS(mfs_head)
        xp, yp = np.array(mfs_wcs.celestial.world_to_pixel(coord)).round().astype(int)
        tt1_p = mfs_i_1[yp, xp]
        tt0_p = mfs_i_0[yp, xp]

        alpha = -1 * tt1_p / tt0_p
        amplitude = tt0_p
        x_0 = mfs_head["RESTFREQ"]

        logger.debug(f"alpha is {alpha}")
        model_I = models.PowerLaw1D(amplitude=amplitude, x_0=x_0, alpha=alpha)

        return StokesIFitResult(
            alpha=alpha,
            amplitude=amplitude,
            x_0=x_0,
            model_repr=model_I.__repr__(),
            modStokesI=model_I(freq),
            fit_dict=None,
        )

    elif do_own_fit:
        logger.info(f"Doing own fit")
        fit_dict = fit_pl(freq=freq, flux=iarr, fluxerr=rmsi, nterms=abs(polyOrd))

        return StokesIFitResult(
            alpha=None,
            amplitude=None,
            x_0=None,
            model_repr=None,
            modStokesI=fit_dict["best_m"],
            fit_dict=fit_dict,
        )

    else:
        return StokesIFitResult(
            alpha=None,
            amplitude=None,
            x_0=None,
            model_repr=None,
            modStokesI=None,
            fit_dict=None,
        )


def update_rmtools_dict(
    mDict: dict,
    fit_dict: dict,
) -> dict:
    """Update the RM-Tools dictionary with the fit results from the Stokes I fit

    Args:
        mDict (dict): The RM-Tools dictionary
        fit_dict (dict): The fit results dictionary

    Returns:
        dict: The updated RM-Tools dictionary
    """
    # Wrangle into format that matches RM-Tools
    mDict["polyCoeffs"] = ",".join(
        [
            # Pad with zeros to length 5
            str(i)
            for i in np.pad(
                fit_dict["best_p"],
                (0, 5 - len(fit_dict["best_p"])),
                "constant",
                constant_values=np.nan,
            )[::-1]
        ]
    )
    mDict["polyCoefferr"] = ",".join(
        [
            str(i)
            for i in np.pad(
                fit_dict["best_e"],
                (0, 5 - len(fit_dict["best_e"])),
                "constant",
                constant_values=np.nan,
            )[::-1]
        ]
    )
    mDict["polyOrd"] = (
        int(fit_dict["best_n"]) if np.isfinite(fit_dict["best_n"]) else float(np.nan)
    )
    mDict["poly_reffreq"] = float(fit_dict["ref_nu"])
    mDict["IfitChiSqRed"] = float(fit_dict["chi_sq_red"])
    for key, val in fit_dict["fit_flag"].items():
        mDict[f"fit_flag_{key}"] = val


@delayed
def rmsynthoncut1d(
    comp: dict,
    beam: dict,
    outdir: str,
    freq: np.ndarray,
    field: str,
    polyOrd: int = 3,
    phiMax_radm2: Optional[float] = None,
    dPhi_radm2: Optional[float] = None,
    nSamples: int = 5,
    weightType: str = "variance",
    fitRMSF: bool = True,
    noStokesI: bool = False,
    showPlots: bool = False,
    savePlots: bool = False,
    debug: bool = False,
    rm_verbose: bool = False,
    fit_function: str = "log",
    tt0: Optional[str] = None,
    tt1: Optional[str] = None,
    ion: bool = False,
    do_own_fit: bool = False,
) -> pymongo.UpdateOne:
    """1D RM synthesis

    Args:
        comp_id (str): RACS component ID
        outdir (str): Output directory
        freq (list): Frequencies in Hz
        host (str): MongoDB host
        field (str): RACS field
        database (bool, optional): Update MongoDB. Defaults to False.
        polyOrd (int, optional): Order of fit to I. Defaults to 3.
        phiMax_radm2 (float, optional): Max FD. Defaults to None.
        dPhi_radm2 (float, optional): Delta FD. Defaults to None.
        nSamples (int, optional): Samples across RMSF. Defaults to 5.
        weightType (str, optional): Weight type. Defaults to 'variance'.
        fitRMSF (bool, optional): Fit RMSF. Defaults to False.
        noStokesI (bool, optional): Ignore Stokes I. Defaults to False.
        showPlots (bool, optional): Show plots. Defaults to False.
        savePlots (bool, optional): Save plots. Defaults to False.
        debug (bool, optional): Turn on debug plots. Defaults to False.
        rm_verbose (bool, optional): Verbose RMsynth. Defaults to False.
    """
    logger.setLevel(logging.INFO)

    iname = comp["Source_ID"]
    cname = comp["Gaussian_ID"]
    ra = comp["RA"]
    dec = comp["Dec"]
    coord = SkyCoord(ra * u.deg, dec * u.deg)
    field_dict = beam["beams"][field]

    # Extract the spectra from the cubelets
    stokes_spectra = extract_all_spectra(
        coord=coord,
        ion=ion,
        field_dict=field_dict,
        outdir=outdir,
    )
    # Check for totally bad data
    for spectrum in stokes_spectra:
        if np.isnan(spectrum.data).all():
            logger.critical(f"Entire data is NaN for {iname} in {spectrum.filename}")
            myquery = {"Gaussian_ID": cname}
            badvalues = {"$set": {"rmsynth1d": False}}
            return pymongo.UpdateOne(myquery, badvalues)

    prefix = f"{os.path.dirname(stokes_spectra.i.filename)}/{cname}"

    # Filter by RMS for outlier rejection
    filtered_stokes_spectra = sigma_clip_spectra(stokes_spectra)

    stokes_i_fit_result = fit_stokes_I(
        freq=freq,
        coord=coord,
        tt0=tt0,
        tt1=tt1,
        do_own_fit=do_own_fit,
        iarr=filtered_stokes_spectra.i.data,
        rmsi=filtered_stokes_spectra.i.rms,
        polyOrd=polyOrd,
    )

    # Check for totally bad data in Q and U
    if (
        np.sum(np.isfinite(filtered_stokes_spectra.q.data)) < 2
        or np.sum(np.isfinite(filtered_stokes_spectra.u.data)) < 2
    ):
        logger.critical(f"{cname} QU data is all NaNs.")
        myquery = {"Gaussian_ID": cname}
        badvalues = {"$set": {"rmsynth1d": False}}
        return pymongo.UpdateOne(myquery, badvalues)
    # And I
    if np.isnan(filtered_stokes_spectra.i.data).all():
        logger.critical(f"{cname} I data is all NaNs.")
        myquery = {"Gaussian_ID": cname}
        badvalues = {"$set": {"rmsynth1d": False}}
        return pymongo.UpdateOne(myquery, badvalues)

    data = [np.array(freq)]
    for stokes in "iqu":
        if noStokesI and stokes == "i":
            continue
        data.append(filtered_stokes_spectra.__getattribute__(stokes).data)
        data.append(filtered_stokes_spectra.__getattribute__(stokes).rms)

    # Run 1D RM-synthesis on the spectra
    np.savetxt(f"{prefix}.dat", np.vstack(data).T, delimiter=" ")
    try:
        logger.info(f"Using {fit_function} to fit Stokes I")
        mDict, aDict = do_RMsynth_1D.run_rmsynth(
            data=data,
            polyOrd=polyOrd,
            phiMax_radm2=phiMax_radm2,
            dPhi_radm2=dPhi_radm2,
            nSamples=nSamples,
            weightType=weightType,
            fitRMSF=fitRMSF,
            noStokesI=noStokesI,
            modStokesI=stokes_i_fit_result.modStokesI,
            nBits=32,
            saveFigures=savePlots,
            showPlots=showPlots,
            verbose=rm_verbose,
            debug=debug,
            fit_function=fit_function,
            prefixOut=prefix,
        )
    except Exception as err:
        traceback.print_tb(err.__traceback__)
        raise err
    if savePlots:
        plt.close("all")
        plotdir = os.path.join(outdir, "plots")
        plot_files = glob(
            os.path.join(os.path.dirname(filtered_stokes_spectra.i.filename), "*.pdf")
        )
        for src in plot_files:
            base = os.path.basename(src)
            dst = os.path.join(plotdir, base)
            copyfile(src, dst)

    # Update model values if own fit was used
    if do_own_fit:
        mDict = update_rmtools_dict(mDict, stokes_i_fit_result.fit_dict)
    else:
        # 0: Improper input parameters (not sure what would trigger this in RM-Tools?)
        # 1-4: One or more of the convergence criteria was met.
        # 5: Reached maximum number of iterations before converging.
        # 6-8: User defined limits for convergence are too small (should not occur, since RM-Tools uses default values)
        # 9: fit failed, reason unknown
        # 16: a fit parameter has become infinite/numerical overflow
        # +64 (can be added to other flags): model gives Stokes I values with S:N < 1 for at least one channel
        # +128 (can be added to other flags): model gives Stokes I values < 0 for at least one channel
        mDict["fit_flag_is_negative"] = mDict["IfitStat"] >= 128
        mDict["fit_flag_is_close_to_zero"] = mDict["IfitStat"] >= 64
        mDict["fit_flag_is_not_finite"] = mDict["IfitStat"] >= 16
        mDict["fit_flag_is_not_normal"] = mDict["IfitStat"] >= 5

    # Ensure JSON serializable
    for k, v in mDict.items():
        if isinstance(v, np.float_):
            mDict[k] = float(v)
        elif isinstance(v, np.float32):
            mDict[k] = float(v)
        elif isinstance(v, np.int_):
            mDict[k] = int(v)
        elif isinstance(v, np.int32):
            mDict[k] = int(v)
        elif isinstance(v, np.ndarray):
            mDict[k] = v.tolist()
        elif isinstance(v, np.bool_):
            mDict[k] = bool(v)

    do_RMsynth_1D.saveOutput(mDict, aDict, prefix, rm_verbose)

    myquery = {"Gaussian_ID": cname}

    # Prep header
    head_dict = dict(filtered_stokes_spectra.i.header)
    head_dict.pop("", None)
    if "COMMENT" in head_dict.keys():
        head_dict["COMMENT"] = str(head_dict["COMMENT"])

    outer_dir = os.path.basename(os.path.dirname(filtered_stokes_spectra.i.filename))

    # Fix for json encoding

    newvalues = {
        "$set": {
            "rm1dfiles": {
                "FDF_dirty": os.path.join(outer_dir, f"{cname}_FDFdirty.dat"),
                "RMSF": os.path.join(outer_dir, f"{cname}_RMSF.dat"),
                "weights": os.path.join(outer_dir, f"{cname}_weight.dat"),
                "summary_dat": os.path.join(outer_dir, f"{cname}_RMsynth.dat"),
                "summary_json": os.path.join(outer_dir, f"{cname}_RMsynth.json"),
            },
            "rmsynth1d": True,
            "header": head_dict,
            "rmsynth_summary": mDict,
            "spectra": {
                "freq": np.array(freq).tolist(),
                "I_model": stokes_i_fit_result.modStokesI.tolist()
                if stokes_i_fit_result.modStokesI is not None
                else None,
                "I_model_params": {
                    "alpha": float(stokes_i_fit_result.alpha)
                    if stokes_i_fit_result.alpha is not None
                    else None,
                    "amplitude": float(stokes_i_fit_result.amplitude)
                    if stokes_i_fit_result.amplitude is not None
                    else None,
                    "x_0": float(stokes_i_fit_result.x_0)
                    if stokes_i_fit_result.x_0 is not None
                    else None,
                    "model_repr": stokes_i_fit_result.model_repr,
                },
                "I": filtered_stokes_spectra.i.data.tolist(),
                "Q": filtered_stokes_spectra.q.data.tolist(),
                "U": filtered_stokes_spectra.u.data.tolist(),
                "I_err": filtered_stokes_spectra.i.rms.tolist(),
                "Q_err": filtered_stokes_spectra.q.rms.tolist(),
                "U_err": filtered_stokes_spectra.u.rms.tolist(),
            },
        }
    }
    return pymongo.UpdateOne(myquery, newvalues)


@delayed
def rmsynthoncut_i(
    comp_id: str,
    outdir: str,
    freq: np.ndarray,
    host: str,
    field: str,
    epoch: int,
    username: Union[str, None] = None,
    password: Union[str, None] = None,
    nSamples: int = 5,
    phiMax_radm2: Union[float, None] = None,
    verbose: bool = False,
    rm_verbose: bool = False,
):
    """RMsynth on Stokes I

    Args:
        comp_id (str): RACS component ID
        freq (list): Frequencies in Hz
        host (str): MongoDB host
        field (str): RACS field
        nSamples ([type]): Samples across the RMSF
        phiMax_radm2 (float): Max FD
        verbose (bool, optional): Verbose output Defaults to False.
        rm_verbose (bool, optional): Verbose RMsynth. Defaults to False.
    """
    logger.setLevel(logging.INFO)

    beams_col, island_col, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )

    # Basic querey
    myquery = {"Gaussian_ID": comp_id}
    doc = comp_col.find_one(myquery)

    if doc is None:
        raise ValueError(f"Component {comp_id} not found")

    iname = doc["Source_ID"]
    cname = doc["Gaussian_ID"]

    beams = beams_col.find_one({"Source_ID": iname})
    if beams is None:
        raise ValueError(f"Beams for {iname} not found")

    ifile = os.path.join(outdir, beams["beams"][field]["i_file"])
    outdir = os.path.dirname(ifile)

    header, dataI = do_RMsynth_3D.readFitsCube(ifile, rm_verbose)

    prefix = f"{outdir}/validation_{cname}"
    # Get source peak from Selavy
    ra = doc["RA"]
    dec = doc["Dec"]
    if len(dataI.shape) == 4:
        # drop Stokes axis
        wcs = WCS(header).dropaxis(2)
    else:
        wcs = WCS(header)

    x, y, z = (
        np.array(wcs.all_world2pix(ra, dec, np.nanmean(freq), 0)).round().astype(int)
    )

    mom = np.nansum(dataI, axis=0)

    plt.ion()
    plt.figure()
    plt.imshow(mom, origin="lower", cmap="cubehelix_r")
    plt.scatter(x, y, c="r", marker="x")
    plt.show()
    _ = input("Press [enter] to continue")  # wait for input from the user
    plt.close()  # close the figure to show the next one.

    data = np.nansum(dataI[:, y - 1 : y + 1 + 1, x - 1 : x + 1 + 1], axis=(1, 2))

    bkgi, rmsi = cubelet_bane(dataI, header)
    rmsi[rmsi == 0] = np.nan
    rmsi[np.isnan(rmsi)] = np.nanmedian(rmsi)
    noise = rmsi

    plt.ion()
    plt.figure()
    plt.step(freq / 1e9, data)

    imod, qArr, uArr, dqArr, duArr, fitDict = create_frac_spectra(
        freqArr=freq,
        IArr=data,
        QArr=data,
        UArr=data,
        dIArr=noise,
        dQArr=noise,
        dUArr=noise,
        polyOrd=3,
        verbose=True,
        debug=False,
    )
    plt.plot(freq / 1e9, imod)
    plt.xlabel(r"$\nu$ [GHz]")
    plt.ylabel(r"$I$ [Jy/beam]")
    plt.tight_layout()
    plt.show()
    _ = input("Press [enter] to continue")  # wait for input from the user
    plt.close()  # close the figure to show the next one.

    data = data - imod
    data = data - np.nanmean(data)
    plt.ion()
    plt.figure()
    plt.step(freq / 1e9, data)
    plt.xlabel(r"$\nu$ [GHz]")
    plt.ylabel(r"$I-\mathrm{model}(I)-\mathrm{mean}(\mathrm{model}(I))$")
    plt.tight_layout()
    plt.show()
    _ = input("Press [enter] to continue")  # wait for input from the user
    plt.close()  # close the figure to show the next one.

    datalist = [freq, data, data, dqArr, duArr]

    phi_max = phiMax_radm2
    mDict, aDict = do_RMsynth_1D.run_rmsynth(
        datalist, phiMax_radm2=phi_max, nSamples=nSamples, verbose=True
    )
    plt.ion()
    plt.figure()
    plt.plot(aDict["phiArr_radm2"], abs(aDict["dirtyFDF"]))
    plt.xlabel(r"$\phi$ [rad m$^{-2}$]")
    plt.ylabel(r"Dirty FDF (Stokes I)")
    plt.tight_layout()
    plt.show()
    _ = input("Press [enter] to continue")  # wait for input from the user
    plt.close()  # close the figure to show the next one.
    do_RMsynth_1D.saveOutput(mDict, aDict, prefix, verbose=verbose)


def main(
    field: str,
    outdir: Path,
    host: str,
    epoch: int,
    username: Union[str, None] = None,
    password: Union[str, None] = None,
    dimension: str = "1d",
    verbose: bool = True,
    database: bool = False,
    validate: bool = False,
    limit: Union[int, None] = None,
    savePlots: bool = False,
    weightType: str = "variance",
    fitRMSF: bool = True,
    phiMax_radm2: Union[float, None] = None,
    dPhi_radm2: Union[float, None] = None,
    nSamples: int = 5,
    polyOrd: int = 3,
    noStokesI: bool = False,
    showPlots: bool = False,
    not_RMSF: bool = False,
    rm_verbose: bool = False,
    debug: bool = False,
    fit_function: str = "log",
    tt0: Union[str, None] = None,
    tt1: Union[str, None] = None,
    ion: bool = False,
    do_own_fit: bool = False,
) -> None:
    outdir = os.path.abspath(outdir)
    outdir = os.path.join(outdir, "cutouts")

    if savePlots:
        plotdir = os.path.join(outdir, "plots")
        try_mkdir(plotdir)

    beams_col, island_col, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )

    beam_query = {"$and": [{f"beams.{field}": {"$exists": True}}]}

    beams = pd.DataFrame(list(beams_col.find(beam_query).sort("Source_ID")))
    beams.set_index("Source_ID", drop=False, inplace=True)
    island_ids = sorted(beams_col.distinct("Source_ID", beam_query))

    isl_query = {"Source_ID": {"$in": island_ids}}
    components = pd.DataFrame(
        list(
            comp_col.find(
                isl_query,
                # Only get required values
                {
                    "Source_ID": 1,
                    "Gaussian_ID": 1,
                    "RA": 1,
                    "Dec": 1,
                },
            ).sort("Source_ID")
        )
    )
    components.set_index("Source_ID", drop=False, inplace=True)
    component_ids = list(components["Gaussian_ID"])

    n_comp = comp_col.count_documents(isl_query)
    n_island = island_col.count_documents(isl_query)

    # Unset rmsynth in db
    if dimension == "1d":
        query_1d = {"$and": [{"Source_ID": {"$in": island_ids}}, {"rmsynth1d": True}]}

        comp_col.update_many(query_1d, {"$set": {"rmsynth1d": False}})

    elif dimension == "3d":
        query_3d = {"$and": [{"Source_ID": {"$in": island_ids}}, {"rmsynth3d": True}]}

        island_col.update(query_3d, {"$set": {"rmsynth3d": False}})

    if limit is not None:
        n_comp = limit
        n_island = limit
        island_ids = island_ids[:limit]
        component_ids = component_ids[:limit]

    # Make frequency file
    freq, freqfile = getfreq(
        os.path.join(outdir, f"{beams.iloc[0]['beams'][f'{field}']['q_file']}"),
        outdir=outdir,
        filename="frequencies.txt",
    )
    freq = np.array(freq)

    outputs = []

    if validate:
        logger.info(f"Running RMsynth on {n_comp} components")
        # We don't run this in parallel!
        for i, comp_id in enumerate(component_ids):
            output = rmsynthoncut_i(
                comp_id=comp_id,
                outdir=outdir,
                freq=freq,
                host=host,
                epoch=epoch,
                field=field,
                username=username,
                password=password,
                nSamples=nSamples,
                phiMax_radm2=phiMax_radm2,
                verbose=verbose,
                rm_verbose=rm_verbose,
            )
            output.compute()

    elif dimension == "1d":
        logger.info(f"Running RMsynth on {n_comp} components")
        for i, (_, comp) in tqdm(
            enumerate(components.iterrows()),
            total=n_comp,
            disable=(not verbose),
            desc="Constructing 1D RMsynth jobs",
        ):
            if i > n_comp + 1:
                break
            else:
                beam = dict(beams.loc[comp["Source_ID"]])
                output = rmsynthoncut1d(
                    comp=comp,
                    beam=beam,
                    outdir=outdir,
                    freq=freq,
                    field=field,
                    polyOrd=polyOrd,
                    phiMax_radm2=phiMax_radm2,
                    dPhi_radm2=dPhi_radm2,
                    nSamples=nSamples,
                    weightType=weightType,
                    fitRMSF=fitRMSF,
                    noStokesI=noStokesI,
                    showPlots=showPlots,
                    savePlots=savePlots,
                    debug=debug,
                    rm_verbose=rm_verbose,
                    fit_function=fit_function,
                    tt0=tt0,
                    tt1=tt1,
                    ion=ion,
                    do_own_fit=do_own_fit,
                )
                outputs.append(output)

    elif dimension == "3d":
        logger.info(f"Running RMsynth on {n_island} islands")

        for i, island_id in enumerate(island_ids):
            if i > n_island + 1:
                break
            else:
                beam = dict(beams.loc[island_id])
                output = rmsynthoncut3d(
                    island_id=island_id,
                    beam=beam,
                    outdir=outdir,
                    freq=freq,
                    field=field,
                    phiMax_radm2=phiMax_radm2,
                    dPhi_radm2=dPhi_radm2,
                    nSamples=nSamples,
                    weightType=weightType,
                    fitRMSF=fitRMSF,
                    not_RMSF=not_RMSF,
                    rm_verbose=rm_verbose,
                    ion=ion,
                )
                outputs.append(output)
    else:
        raise ValueError("An incorrect RMSynth mode has been configured. ")

    futures = chunk_dask(
        outputs=outputs,
        task_name="RMsynth",
        progress_text="Running RMsynth",
        verbose=verbose,
    )

    if database:
        logger.info("Updating database...")
        updates = [f.compute() for f in futures]
        # Remove None values
        updates = [u for u in updates if u is not None]
        logger.info("Sending updates to database...")
        if dimension == "1d":
            db_res = comp_col.bulk_write(updates, ordered=False)
            logger.info(pformat(db_res.bulk_api_result))
        elif dimension == "3d":
            db_res = island_col.bulk_write(updates, ordered=False)
            logger.info(pformat(db_res.bulk_api_result))
    logger.info("RMsynth done!")


def cli():
    """Command-line interface"""
    import argparse

    from astropy.utils.exceptions import AstropyWarning

    warnings.simplefilter("ignore", category=AstropyWarning)
    from astropy.io.fits.verify import VerifyWarning

    warnings.simplefilter("ignore", category=VerifyWarning)
    warnings.simplefilter("ignore", category=RuntimeWarning)
    # Help string to be shown using the -h option
    descStr = f"""
    {logo_str}
    Arrakis Stage 5:
    Run RMsynthesis on cubelets.

    Note: Runs on brightest sources first.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field to mosaic - e.g. 2132-50A."
    )
    parser.add_argument(
        "outdir",
        metavar="outdir",
        type=Path,
        help="Directory containing cutouts (in subdir outdir/cutouts).",
    )

    parser.add_argument(
        "host",
        metavar="host",
        type=str,
        help="Host of mongodb (probably $hostname -i).",
    )

    parser.add_argument(
        "-e",
        "--epoch",
        type=int,
        default=0,
        help="Epoch of observation.",
    )

    parser.add_argument(
        "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )

    parser.add_argument(
        "--dimension",
        dest="dimension",
        default="1d",
        help="How many dimensions for RMsynth [1d] or '3d'.",
    )

    parser.add_argument(
        "-v", dest="verbose", action="store_true", help="verbose output [False]."
    )

    parser.add_argument(
        "--ion", action="store_true", help="Use ionospheric-corrected data [False]."
    )

    parser.add_argument(
        "-m", dest="database", action="store_true", help="Add data to MongoDB [False]."
    )

    parser.add_argument(
        "--tt0",
        default=None,
        type=str,
        help="TT0 MFS image -- will be used for model of Stokes I -- also needs --tt1.",
    )

    parser.add_argument(
        "--tt1",
        default=None,
        type=str,
        help="TT1 MFS image -- will be used for model of Stokes I -- also needs --tt0.",
    )

    parser.add_argument(
        "--validate",
        dest="validate",
        action="store_true",
        help="Run on Stokes I [False].",
    )

    parser.add_argument(
        "--limit",
        dest="limit",
        default=None,
        type=int,
        help="Limit number of sources [All].",
    )
    parser.add_argument(
        "--own_fit",
        dest="do_own_fit",
        action="store_true",
        help="Use own Stokes I fit function [False].",
    )

    # RM-tools args
    parser.add_argument(
        "-sp", "--savePlots", action="store_true", help="save the plots [False]."
    )
    parser.add_argument(
        "-w",
        dest="weightType",
        default="variance",
        help="weighting [variance] (all 1s) or 'uniform'.",
    )
    parser.add_argument(
        "-f",
        dest="fit_function",
        type=str,
        default="log",
        help="Stokes I fitting function: 'linear' or ['log'] polynomials.",
    )
    parser.add_argument(
        "-t",
        dest="fitRMSF",
        action="store_true",
        help="Fit a Gaussian to the RMSF [False]",
    )
    parser.add_argument(
        "-l",
        dest="phiMax_radm2",
        type=float,
        default=None,
        help="Absolute max Faraday depth sampled (overrides NSAMPLES) [Auto].",
    )
    parser.add_argument(
        "-d",
        dest="dPhi_radm2",
        type=float,
        default=None,
        help="Width of Faraday depth channel [Auto].",
    )
    parser.add_argument(
        "-s",
        dest="nSamples",
        type=float,
        default=5,
        help="Number of samples across the FWHM RMSF.",
    )
    parser.add_argument(
        "-o",
        dest="polyOrd",
        type=int,
        default=3,
        help="polynomial order to fit to I spectrum [3].",
    )
    parser.add_argument(
        "-i",
        dest="noStokesI",
        action="store_true",
        help="ignore the Stokes I spectrum [False].",
    )
    parser.add_argument(
        "-p", dest="showPlots", action="store_true", help="show the plots [False]."
    )
    parser.add_argument(
        "-R",
        dest="not_RMSF",
        action="store_true",
        help="Skip calculation of RMSF? [False]",
    )
    parser.add_argument(
        "-rmv", dest="rm_verbose", action="store_true", help="Verbose RMsynth [False]."
    )
    parser.add_argument(
        "-D",
        dest="debug",
        action="store_true",
        help="turn on debugging messages & plots [False].",
    )

    args = parser.parse_args()

    if args.tt0 and not args.tt1:
        parser.error("the following arguments are required: tt1")
    elif args.tt1 and not args.tt0:
        parser.error("the following arguments are required: tt0")

    verbose = args.verbose
    rmv = args.rm_verbose
    if rmv:
        logger.setLevel(logger.DEBUG)
    elif verbose:
        logger.setLevel(logger.INFO)

    cluster = LocalCluster(
        # n_workers=12, processes=True, threads_per_worker=1,
        local_directory="/dev/shm"
    )
    client = Client(cluster)
    logger.debug(client)

    test_db(
        host=args.host, username=args.username, password=args.password, verbose=verbose
    )

    main(
        field=args.field,
        outdir=Path(args.outdir),
        host=args.host,
        username=args.username,
        password=args.password,
        dimension=args.dimension,
        verbose=verbose,
        database=args.database,
        validate=args.validate,
        limit=args.limit,
        savePlots=args.savePlots,
        weightType=args.weightType,
        fitRMSF=args.fitRMSF,
        phiMax_radm2=args.phiMax_radm2,
        dPhi_radm2=args.dPhi_radm2,
        nSamples=args.nSamples,
        polyOrd=args.polyOrd,
        noStokesI=args.noStokesI,
        showPlots=args.showPlots,
        not_RMSF=args.not_RMSF,
        rm_verbose=args.rm_verbose,
        debug=args.debug,
        fit_function=args.fit_function,
        tt0=args.tt0,
        tt1=args.tt1,
        ion=args.ion,
        do_own_fit=args.do_own_fit,
    )


if __name__ == "__main__":
    cli()
