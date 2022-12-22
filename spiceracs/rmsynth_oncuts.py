#!/usr/bin/env python3
"""Run RM-CLEAN on cutouts in parallel"""
import functools
import json
import logging as log
import os
import pdb
import subprocess
import sys
import time
import traceback
import warnings
from glob import glob
from pprint import pformat, pprint
from shutil import copyfile
from typing import List, Optional, Tuple, Union

import astropy.units as u
import dask
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
import pymongo
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling import fitting, models
from astropy.stats import mad_std, sigma_clip
from astropy.wcs import WCS
from dask import delayed
from dask.diagnostics import ProgressBar
from dask.distributed import Client, LocalCluster, progress, wait
from IPython import embed
from RMtools_1D import do_RMsynth_1D
from RMtools_3D import do_RMsynth_3D
from RMutils.util_misc import create_frac_spectra
from RMutils.util_plotTk import plot_rmsf_fdf_fig
from spectral_cube import SpectralCube
from tqdm import tqdm, trange

from spiceracs.utils import (
    MyEncoder,
    chunk_dask,
    fit_pl,
    get_db,
    getfreq,
    test_db,
    tqdm_dask,
    try_mkdir,
)


@delayed
def rmsynthoncut3d(
    island_id: str,
    beam: dict,
    outdir: str,
    freq: np.ndarray,
    field: str,
    phiMax_radm2: Union[None, float] = None,
    dPhi_radm2: Union[None, float] = None,
    nSamples: int = 5,
    weightType: str="variance",
    fitRMSF: bool=True,
    not_RMSF: bool=False,
    rm_verbose:bool=False,
    ion: bool=False,
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
        log.critical(f"Cubelet {iname} is entirely NaN")
        myquery = {"Source_ID": iname}
        badvalues = {
            "$set": {
                "rmsynth3d": False,
            }
        }
        return pymongo.UpdateOne(myquery, badvalues)
    rmsi = estimate_noise_annulus(dataI.shape[2] // 2, dataI.shape[1] // 2, dataI)
    rmsi[rmsi == 0] = np.nan
    rmsi[np.isnan(rmsi)] = np.nanmedian(rmsi)

    # rmsq = rms_1d(dataQ)
    rmsq = estimate_noise_annulus(dataQ.shape[2] // 2, dataQ.shape[1] // 2, dataQ)
    rmsq[rmsq == 0] = np.nan
    rmsq[np.isnan(rmsq)] = np.nanmedian(rmsq)

    # rmsu = rms_1d(dataU)
    rmsu = estimate_noise_annulus(dataU.shape[2] // 2, dataU.shape[1] // 2, dataU)
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


@delayed
def rms_1d(data):
    """Compute RMS from bounding pixels"""
    Nfreq, Ndec, Nra = data.shape
    mask = np.ones((Ndec, Nra), dtype=np.bool)
    mask[3:-3, 3:-3] = False
    rms = np.nanstd(data[:, mask], axis=1)
    return rms


def estimate_noise_annulus(x_center, y_center, cube):
    """
    Noise estimation for annulus taken around point source. Annulus has fixed
    inner radius of 10 and outer radius of 31. Function makes an annulus shaped
    mask, then for each source applies the mask at each frequency and takes the
    standard deviation.

    ​Inputs: Array of sets of pixel coordinates (y-position,x-position) for
    sources, Stokes cube (assumes 4 axes), array of flagged channels (can be an
    empty array), number of frequency channels.

    ​Output: 2D array of standard deviation values with shape (length of
    coordinate array, number of unflagged frequency channels).
    """
    cube = np.nan_to_num(cube, nan=0)
    inner_radius = 10
    # Set outer radius to cutout edge if default value is too big
    if min(cube.shape[-2:]) <= 62:
        outer_radius = min(cube.shape[-2:]) // 2 - 1
    else:
        outer_radius = 31

    lenfreq = cube.shape[0]
    naxis = len(cube.shape)
    err = np.zeros(lenfreq)
    # try:
    y, x = np.ogrid[
        -1 * outer_radius : outer_radius + 1, -1 * outer_radius : outer_radius + 1
    ]
    grid_mask = np.logical_or(
        x**2 + y**2 < inner_radius**2, x**2 + y**2 > outer_radius**2
    )
    for i in range(lenfreq):
        if naxis == 4:
            grid = cube[
                i,
                0,
                y_center - outer_radius : y_center + outer_radius + 1,
                x_center - outer_radius : x_center + outer_radius + 1,
            ]
        else:  # naxis ==3
            grid = cube[
                i,
                y_center - outer_radius : y_center + outer_radius + 1,
                x_center - outer_radius : x_center + outer_radius + 1,
            ]

        # Calculate the MADFM, and convert to standard sigma:
        noisepix = np.ma.masked_array(grid, grid_mask)
        # if (noisepix == np.nan).any():
        #    embed
        err[i] = np.ma.median(np.ma.fabs(noisepix - np.ma.median(noisepix))) / 0.6745
    err[err == 0] = np.nan
    return err


@delayed
def rmsynthoncut1d(
    comp: dict,
    beam: dict,
    outdir: str,
    freq: np.ndarray,
    field: str,
    polyOrd: int = 3,
    phiMax_radm2: Union[float,None] = None,
    dPhi_radm2: Union[float,None] = None,
    nSamples: int = 5,
    weightType: str = "variance",
    fitRMSF: bool = True,
    noStokesI: bool = False,
    showPlots: bool = False,
    savePlots: bool = False,
    debug: bool = False,
    rm_verbose: bool = False,
    fit_function: str = "log",
    tt0: Union[str,None] = None,
    tt1: Union[str,None] = None,
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
    iname = comp["Source_ID"]
    cname = comp["Gaussian_ID"]
    ifile = os.path.join(outdir, beam["beams"][field]["i_file"])
    if ion:
        qfile = os.path.join(outdir, beam["beams"][field]["q_file_ion"])
        ufile = os.path.join(outdir, beam["beams"][field]["u_file_ion"])
    else:
        qfile = os.path.join(outdir, beam["beams"][field]["q_file"])
        ufile = os.path.join(outdir, beam["beams"][field]["u_file"])

    header, dataQ = do_RMsynth_3D.readFitsCube(qfile, rm_verbose)
    header, dataU = do_RMsynth_3D.readFitsCube(ufile, rm_verbose)
    header, dataI = do_RMsynth_3D.readFitsCube(ifile, rm_verbose)

    dataQ = np.squeeze(dataQ)
    dataU = np.squeeze(dataU)
    dataI = np.squeeze(dataI)

    if np.isnan(dataI).all() or np.isnan(dataQ).all() or np.isnan(dataU).all():
        log.critical(f"Entire data is NaN for {iname}")
        myquery = {"Gaussian_ID": cname}
        badvalues = {"$set": {"rmsynth1d": False}}
        return pymongo.UpdateOne(myquery, badvalues)

    rmsi = estimate_noise_annulus(dataI.shape[2] // 2, dataI.shape[1] // 2, dataI)
    rmsi[rmsi == 0] = np.nan
    rmsi[np.isnan(rmsi)] = np.nanmedian(rmsi)

    rmsq = estimate_noise_annulus(dataQ.shape[2] // 2, dataQ.shape[1] // 2, dataQ)
    rmsq[rmsq == 0] = np.nan
    rmsq[np.isnan(rmsq)] = np.nanmedian(rmsq)

    rmsu = estimate_noise_annulus(dataU.shape[2] // 2, dataU.shape[1] // 2, dataU)
    rmsu[rmsu == 0] = np.nan
    rmsu[np.isnan(rmsu)] = np.nanmedian(rmsu)

    prefix = f"{os.path.dirname(ifile)}/{cname}"

    ra = comp["RA"]
    dec = comp["Dec"]
    coord = SkyCoord(ra * u.deg, dec * u.deg)
    if len(dataI.shape) == 4:
        # drop Stokes axis
        wcs = WCS(header).dropaxis(2)
    else:
        wcs = WCS(header)

    x, y = np.array(wcs.celestial.world_to_pixel(coord)).round().astype(int)

    qarr = dataQ[:, y, x]
    uarr = dataU[:, y, x]
    iarr = dataI[:, y, x]

    iarr[iarr == 0] = np.nan
    qarr[qarr == 0] = np.nan
    uarr[uarr == 0] = np.nan

    i_filter = sigma_clip(rmsi, sigma=5, stdfunc=mad_std)
    q_filter = sigma_clip(rmsq, sigma=5, stdfunc=mad_std)
    u_filter = sigma_clip(rmsu, sigma=5, stdfunc=mad_std)

    filter_idx = (i_filter.mask) | (q_filter.mask) | (u_filter.mask)

    iarr[filter_idx] = np.nan
    qarr[filter_idx] = np.nan
    uarr[filter_idx] = np.nan
    rmsi[filter_idx] = np.nan
    rmsq[filter_idx] = np.nan
    rmsu[filter_idx] = np.nan

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

        log.debug(f"alpha is {alpha}")
        model_I = models.PowerLaw1D(amplitude=amplitude, x_0=x_0, alpha=alpha)
        modStokesI = model_I(freq)
        model_repr = model_I.__repr__()

    elif do_own_fit:
        log.debug(f"Doing own fit")
        fit_dict = fit_pl(freq=freq, flux=iarr, fluxerr=rmsi, nterms=abs(polyOrd))
        alpha = None
        amplitude = None
        x_0 = None
        model_repr = None
        modStokesI = fit_dict["best_m"]

    else:
        alpha = None
        amplitude = None
        x_0 = None
        model_repr = None
        modStokesI = None

    if np.sum(np.isfinite(qarr)) < 2 or np.sum(np.isfinite(uarr)) < 2:
        log.critical(f"{cname} QU data is all NaNs.")
        myquery = {"Gaussian_ID": cname}
        badvalues = {"$set": {"rmsynth1d": False}}
        return pymongo.UpdateOne(myquery, badvalues)
    if noStokesI:
        data = [np.array(freq), qarr, uarr, rmsq, rmsu]
    else:
        data = [np.array(freq), iarr, qarr, uarr, rmsi, rmsq, rmsu]

    if np.isnan(iarr).all():
        log.critical(f"{cname} I data is all NaNs.")
        myquery = {"Gaussian_ID": cname}
        badvalues = {"$set": {"rmsynth1d": False}}
        return pymongo.UpdateOne(myquery, badvalues)

    # Run 1D RM-synthesis on the spectra
    np.savetxt(f"{prefix}.dat", np.vstack(data).T, delimiter=" ")
    try:
        log.debug(f"Using {fit_function} to fit Stokes I")
        mDict, aDict = do_RMsynth_1D.run_rmsynth(
            data=data,
            polyOrd=polyOrd,
            phiMax_radm2=phiMax_radm2,
            dPhi_radm2=dPhi_radm2,
            nSamples=nSamples,
            weightType=weightType,
            fitRMSF=fitRMSF,
            noStokesI=noStokesI,
            modStokesI=modStokesI,
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
        plot_files = glob(os.path.join(os.path.dirname(ifile), "*.pdf"))
        for src in plot_files:
            base = os.path.basename(src)
            dst = os.path.join(plotdir, base)
            copyfile(src, dst)

    # Update model values if own fit was used
    if do_own_fit:
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
            int(fit_dict["best_n"])
            if np.isfinite(fit_dict["best_n"])
            else float(np.nan)
        )
        mDict["poly_reffreq"] = float(fit_dict["ref_nu"])
        mDict["IfitChiSqRed"] = float(fit_dict["chi_sq_red"])
        for key, val in fit_dict["fit_flag"].items():
            mDict[f"fit_flag_{key}"] = val
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
    head_dict = dict(header)
    head_dict.pop("", None)
    head_dict["COMMENT"] = str(head_dict["COMMENT"])

    outer_dir = os.path.basename(os.path.dirname(ifile))

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
                "I_model": modStokesI.tolist() if modStokesI is not None else None,
                "I_model_params": {
                    "alpha": float(alpha) if alpha is not None else None,
                    "amplitude": float(amplitude) if amplitude is not None else None,
                    "x_0": float(x_0) if x_0 is not None else None,
                    "model_repr": model_repr,
                },
                "I": iarr.tolist(),
                "Q": qarr.tolist(),
                "U": uarr.tolist(),
                "I_err": rmsi.tolist(),
                "Q_err": rmsq.tolist(),
                "U_err": rmsu.tolist(),
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
    username:Union[str,None]=None,
    password:Union[str,None]=None,
    nSamples:int=5,
    phiMax_radm2:Union[float,None]=None,
    verbose:bool=False,
    rm_verbose:bool=False,
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
    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )

    # Basic querey
    myquery = {"Gaussian_ID": comp_id}
    doc = comp_col.find_one(myquery)

    iname = doc["Source_ID"]
    cname = doc["Gaussian_ID"]

    beams = beams_col.find_one({"Source_ID": iname})
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

    rmsi = estimate_noise_annulus(dataI.shape[2] // 2, dataI.shape[1] // 2, dataI)
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
    outdir: str,
    host: str,
    client: Client,
    username: Union[str,None] = None,
    password: Union[str,None] = None,
    dimension: str = "1d",
    verbose: bool = True,
    database: bool = False,
    validate: bool = False,
    limit: Union[int,None] = None,
    savePlots: bool = False,
    weightType: str = "variance",
    fitRMSF: bool = True,
    phiMax_radm2: Union[float,None] = None,
    dPhi_radm2: Union[float,None] = None,
    nSamples: int = 5,
    polyOrd: int = 3,
    noStokesI: bool = False,
    showPlots: bool = False,
    not_RMSF: bool = False,
    rm_verbose: bool = False,
    debug: bool = False,
    fit_function: str = "log",
    tt0: Union[str,None] = None,
    tt1: Union[str,None] = None,
    ion: bool = False,
    do_own_fit: bool = False,
) -> None:

    outdir = os.path.abspath(outdir)
    outdir = os.path.join(outdir, "cutouts")

    if savePlots:
        plotdir = os.path.join(outdir, "plots")
        try_mkdir(plotdir)

    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )

    beam_query = {
        "$and": [{f"beams.{field}": {"$exists": True}}, {f"beams.{field}.DR1": True}]
    }

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
        log.info(f"Running RMsynth on {n_comp} components")
        # We don't run this in parallel!
        for i, comp_id in enumerate(component_ids):
            output = rmsynthoncut_i(
                comp_id=comp_id,
                outdir=outdir,
                freq=freq,
                host=host,
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
        log.info(f"Running RMsynth on {n_comp} components")
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
        log.info(f"Running RMsynth on {n_island} islands")

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

    futures = chunk_dask(
        outputs=outputs,
        client=client,
        task_name="RMsynth",
        progress_text="Running RMsynth",
        verbose=verbose,
    )

    if database:
        log.info("Updating database...")
        updates = [f.compute() for f in futures]
        # Remove None values
        updates = [u for u in updates if u is not None]
        log.info("Sending updates to database...")
        if dimension == "1d":
            db_res = comp_col.bulk_write(updates, ordered=False)
            log.info(pformat(db_res.bulk_api_result))
        elif dimension == "3d":
            db_res = island_col.bulk_write(updates, ordered=False)
            log.info(pformat(db_res.bulk_api_result))
    log.info("RMsynth done!")


def cli():
    """Command-line interface"""
    import argparse

    from astropy.utils.exceptions import AstropyWarning

    warnings.simplefilter("ignore", category=AstropyWarning)
    from astropy.io.fits.verify import VerifyWarning

    warnings.simplefilter("ignore", category=VerifyWarning)
    warnings.simplefilter("ignore", category=RuntimeWarning)
    # Help string to be shown using the -h option
    logostr = """
     mmm   mmm   mmm   mmm   mmm
     )-(   )-(   )-(   )-(   )-(
    ( S ) ( P ) ( I ) ( C ) ( E )
    |   | |   | |   | |   | |   |
    |___| |___| |___| |___| |___|
     mmm     mmm     mmm     mmm
     )-(     )-(     )-(     )-(
    ( R )   ( A )   ( C )   ( S )
    |   |   |   |   |   |   |   |
    |___|   |___|   |___|   |___|

    """

    # Help string to be shown using the -h option
    descStr = f"""
    {logostr}
    SPICE-RACS Stage 5:
    Run RMsynthesis on cubelets.

    Note: Runs on brightest sources first.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field to mosaic - e.g. 2132-50A."
    )
    parser.add_argument(
        "outdir",
        metavar="outdir",
        type=str,
        help="Directory containing cutouts (in subdir outdir/cutouts).",
    )

    parser.add_argument(
        "host",
        metavar="host",
        type=str,
        help="Host of mongodb (probably $hostname -i).",
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
        log.basicConfig(
            level=log.DEBUG,
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    elif verbose:
        log.basicConfig(
            level=log.INFO,
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    else:
        log.basicConfig(
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

    cluster = LocalCluster(
        # n_workers=12, processes=True, threads_per_worker=1,
        local_directory="/dev/shm"
    )
    client = Client(cluster)
    log.debug(client)

    test_db(
        host=args.host, username=args.username, password=args.password, verbose=verbose
    )

    main(
        field=args.field,
        outdir=args.outdir,
        host=args.host,
        username=args.username,
        password=args.password,
        client=client,
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
