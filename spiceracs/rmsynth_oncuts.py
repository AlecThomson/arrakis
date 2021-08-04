#!/usr/bin/env python3
from spiceracs.utils import getfreq, MyEncoder, tqdm_dask, try_mkdir
import json
import numpy as np
import os
from glob import glob
from shutil import copyfile
import pymongo
import sys
import subprocess
import time
from tqdm import tqdm, trange
import warnings
from RMtools_1D import do_RMsynth_1D
from RMtools_3D import do_RMsynth_3D
from spectral_cube import SpectralCube
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clip, mad_std
from astropy.coordinates import SkyCoord
from astropy.modeling import models, fitting
import astropy.units as u
import matplotlib.pyplot as plt
from RMutils.util_plotTk import plot_rmsf_fdf_fig
from RMutils.util_misc import create_frac_spectra
import functools
import psutil
import pdb
from IPython import embed
import dask
from dask import delayed
from dask.distributed import Client, progress, LocalCluster, wait
from dask.diagnostics import ProgressBar
import traceback



# def safe_mongocall(call):
#     def _safe_mongocall(*args, **kwargs):
#         for i in range(5):
#             try:
#                 return call(*args, **kwargs)
#             except pymongo.errors.AutoReconnect:
#                 time.sleep(np.random.random() / 100)
#         print('Error: Failed operation!')
#     return _safe_mongocall


@delayed
def rmsynthoncut3d(island_id,
                   outdir,
                   freq,
                   host,
                   field,
                   database=False,
                   phiMax_radm2=None,
                   dPhi_radm2=None,
                   nSamples=5,
                   weightType='variance',
                   fitRMSF=True,
                   not_RMSF=False,
                   rm_verbose=False,
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
    # default connection (ie, local)
    with pymongo.MongoClient(host=host, connect=False) as dbclient:
        mydb = dbclient['spiceracs']  # Create/open database
        isl_col = mydb['islands']  # Create/open collection
        comp_col = mydb['components']  # Create/open collection
        beams_col = mydb['beams']  # Create/open collection

    # Basic querey
    myquery = {"Source_ID": island_id}

    doc = comp_col.find_one(myquery)
    iname = doc['Source_ID']
    beam = beams_col.find_one({'Source_ID': iname})

    ifile = os.path.join(outdir, beam['beams'][field]['i_file'])
    qfile = os.path.join(outdir, beam['beams'][field]['q_file'])
    ufile = os.path.join(outdir, beam['beams'][field]['u_file'])
    # vfile = beam['beams'][field]['v_file']

    header, dataQ = do_RMsynth_3D.readFitsCube(qfile, rm_verbose)
    header, dataU = do_RMsynth_3D.readFitsCube(ufile, rm_verbose)
    header, dataI = do_RMsynth_3D.readFitsCube(ifile, rm_verbose)

    dataQ = np.squeeze(dataQ)
    dataU = np.squeeze(dataU)
    dataI = np.squeeze(dataI)

    if np.isnan(dataI).all() or np.isnan(dataQ).all() or np.isnan(dataU).all():
        return
        #rmsi = rms_1d(dataI)
    rmsi = estimate_noise_annulus(
        dataI.shape[2]//2,
        dataI.shape[1]//2,
        dataI
    )
    rmsi[rmsi == 0] = np.nan
    rmsi[np.isnan(rmsi)] = np.nanmedian(rmsi)

    #rmsq = rms_1d(dataQ)
    rmsq = estimate_noise_annulus(
        dataQ.shape[2]//2,
        dataQ.shape[1]//2,
        dataQ
    )
    rmsq[rmsq == 0] = np.nan
    rmsq[np.isnan(rmsq)] = np.nanmedian(rmsq)

    #rmsu = rms_1d(dataU)
    rmsu = estimate_noise_annulus(
        dataU.shape[2]//2,
        dataU.shape[1]//2,
        dataU
    )
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
        not_rmsf=not_RMSF
    )

    prefix = f"{iname}_"
    # Write to files
    do_RMsynth_3D.writefits(dataArr,
                            headtemplate=header,
                            fitRMSF=fitRMSF,
                            prefixOut=prefix,
                            outDir=os.path.dirname(ifile),
                            write_seperate_FDF=True,
                            not_rmsf=not_RMSF,
                            nBits=32,
                            verbose=rm_verbose
                            )

    if database:
        myquery = {"Source_ID": iname}
        # Prep header
        head_dict = dict(header)
        head_dict.pop('', None)
        head_dict['COMMENT'] = str(head_dict['COMMENT'])

        outer_dir = os.path.basename(os.path.dirname(ifile))

        newvalues = {
            "$set":
            {
                "rm3dfiles":
                {
                    "FDF_real_dirty": os.path.join(outer_dir, f"{prefix}FDF_real_dirty.fits"),
                    "FDF_im_dirty": os.path.join(outer_dir, f"{prefix}FDF_im_dirty.fits"),
                    "FDF_tot_dirty": os.path.join(outer_dir, f"{prefix}FDF_tot_dirty.fits"),
                    "RMSF_real": os.path.join(outer_dir, f"{prefix}RMSF_real.fits"),
                    "RMSF_tot": os.path.join(outer_dir, f"{prefix}RMSF_tot.fits"),
                    "RMSF_FWHM": os.path.join(outer_dir, f"{prefix}RMSF_FWHM.fits"),
                },
                "rmsynth3d": True,
                "header": dict(header)
            }
        }

        isl_col.update_one(myquery, newvalues)


@delayed
def rms_1d(data):
    """Compute RMS from bounding pixels
    """
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
    y, x = np.ogrid[-1*outer_radius:outer_radius +
                    1, -1*outer_radius:outer_radius+1]
    grid_mask = np.logical_or(
        x**2+y**2 < inner_radius**2, x**2+y**2 > outer_radius**2)
    for i in range(lenfreq):
        if naxis == 4:
            grid = cube[i, 0, y_center-outer_radius:y_center+outer_radius+1,
                        x_center-outer_radius:x_center+outer_radius+1]
        else:  # naxis ==3
            grid = cube[i, y_center-outer_radius:y_center+outer_radius+1,
                        x_center-outer_radius:x_center+outer_radius+1]

        # Calculate the MADFM, and convert to standard sigma:
        noisepix = np.ma.masked_array(grid, grid_mask)
        # if (noisepix == np.nan).any():
        #    embed
        err[i] = np.ma.median(np.ma.fabs(
            noisepix - np.ma.median(noisepix))) / 0.6745
    err[err == 0] = np.nan
    return err
    # except np.ma.core.MaskError:
    #     print('Too many NaNs in cutout - skipping...')
    #     err *= np.nan
    #     return err


@delayed
def rmsynthoncut1d(comp_id,
                   outdir,
                   freq,
                   host,
                   field,
                   database=False,
                   polyOrd=3,
                   phiMax_radm2=None,
                   dPhi_radm2=None,
                   nSamples=5,
                   weightType='variance',
                   fitRMSF=True,
                   noStokesI=False,
                   showPlots=False,
                   savePlots=False,
                   debug=False,
                   rm_verbose=False,
                   fit_function='log',
                   tt0=None,
                   tt1=None
                   ):
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
    # default connection (ie, local)
    with pymongo.MongoClient(host=host, connect=False) as dbclient:
        mydb = dbclient['spiceracs']  # Create/open database
        isl_col = mydb['islands']  # Create/open collection
        comp_col = mydb['components']  # Create/open collection
        beams_col = mydb['beams']  # Create/open collection

        # Basic querey
        myquery = {"Gaussian_ID": comp_id}

        doc = comp_col.find_one(myquery)
        iname = doc['Source_ID']
        cname = doc['Gaussian_ID']
        beam = beams_col.find_one({'Source_ID': iname})

        ifile = os.path.join(outdir, beam['beams'][field]['i_file'])
        qfile = os.path.join(outdir, beam['beams'][field]['q_file'])
        ufile = os.path.join(outdir, beam['beams'][field]['u_file'])
        # vfile = beam['beams'][field]['v_file']

        header, dataQ = do_RMsynth_3D.readFitsCube(qfile, rm_verbose)
        header, dataU = do_RMsynth_3D.readFitsCube(ufile, rm_verbose)
        header, dataI = do_RMsynth_3D.readFitsCube(ifile, rm_verbose)

        dataQ = np.squeeze(dataQ)
        dataU = np.squeeze(dataU)
        dataI = np.squeeze(dataI)

        if np.isnan(dataI).all() or np.isnan(dataQ).all() or np.isnan(dataU).all():
            return
            #rmsi = rms_1d(dataI)
        rmsi = estimate_noise_annulus(
            dataI.shape[2]//2,
            dataI.shape[1]//2,
            dataI
        )
        rmsi[rmsi == 0] = np.nan
        rmsi[np.isnan(rmsi)] = np.nanmedian(rmsi)

        #rmsq = rms_1d(dataQ)
        rmsq = estimate_noise_annulus(
            dataQ.shape[2]//2,
            dataQ.shape[1]//2,
            dataQ
        )
        rmsq[rmsq == 0] = np.nan
        rmsq[np.isnan(rmsq)] = np.nanmedian(rmsq)

        #rmsu = rms_1d(dataU)
        rmsu = estimate_noise_annulus(
            dataU.shape[2]//2,
            dataU.shape[1]//2,
            dataU
        )
        rmsu[rmsu == 0] = np.nan
        rmsu[np.isnan(rmsu)] = np.nanmedian(rmsu)

        # if np.isnan(rmsi).all() or np.isnan(rmsq).all() or np.isnan(rmsu).all():
        #     print(f'RMS data is all NaNs. Skipping component {cname}...')
        #     return

        prefix = f"{os.path.dirname(ifile)}/{cname}"

        ra = doc['RA']
        dec = doc['Dec']
        coord = SkyCoord(ra*u.deg,dec*u.deg)
        if len(dataI.shape) == 4:
            # drop Stokes axis
            wcs = WCS(header).dropaxis(2)
        else:
            wcs = WCS(header)

        x, y = np.array(
            wcs.celestial.world_to_pixel(coord)
        ).round().astype(int)

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
            alpha = tt1_p / tt0_p
            if rm_verbose:
                print(f'alpha is {alpha}')
            model_I = models.PowerLaw1D(tt0_p, mfs_head['RESTFREQ'], alpha=-alpha)
            modStokesI = model_I(freq)

        else:
            modStokesI=None

        if np.sum(np.isfinite(qarr)) < 2 or np.sum(np.isfinite(uarr)) < 2:
            print(f'QU data is all NaNs. Skipping component {cname}...')
            return
        else:
            if noStokesI:
                idx = np.isnan(qarr) | np.isnan(uarr)
                data = [np.array(freq), qarr,
                        uarr, rmsq, rmsu]
            else:
                if np.isnan(iarr).all():
                    print(f'I data is all NaNs. Skipping component {cname}...')
                    return
                else:
                    idx = np.isnan(qarr) | np.isnan(uarr) | np.isnan(iarr)
                    data = [np.array(freq), iarr, qarr,
                            uarr, rmsi, rmsq, rmsu]

            # Run 1D RM-synthesis on the spectra
            np.savetxt(f"{prefix}.dat", np.vstack(data).T, delimiter=' ')
            try:
                mDict, aDict = do_RMsynth_1D.run_rmsynth(data=data,
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
                plotdir = os.path.join(outdir, 'plots')
                plot_files = glob(os.path.join(os.path.dirname(ifile), '*.pdf'))
                for src in plot_files:
                    base = os.path.basename(src)
                    dst = os.path.join(plotdir,base)
                    copyfile(src, dst)

            do_RMsynth_1D.saveOutput(
                mDict, aDict, prefix, rm_verbose)

            if database:
                myquery = {"Gaussian_ID": cname}

                # Prep header
                head_dict = dict(header)
                head_dict.pop('', None)
                head_dict['COMMENT'] = str(head_dict['COMMENT'])

                outer_dir = os.path.basename(os.path.dirname(ifile))

                newvalues = {
                    "$set": {
                        f"rm1dfiles": {
                            "FDF_dirty": os.path.join(outer_dir, f"{cname}_FDFdirty.dat"),
                            "RMSF": os.path.join(outer_dir, f"{cname}_RMSF.dat"),
                            "weights": os.path.join(outer_dir, f"{cname}_weight.dat"),
                            "summary_dat": os.path.join(outer_dir, f"{cname}_RMsynth.dat"),
                            "summary_json": os.path.join(outer_dir, f"{cname}_RMsynth.json"),
                        },
                        f"rmsynth1d": True,
                        "header": head_dict,
                        f"rmsynth_summary": mDict
                    }
                }
                comp_col.update_one(myquery, newvalues)


@delayed
def rmsynthoncut_i(comp_id,
                   outdir,
                   freq,
                   host,
                   field,
                   nSamples=5,
                   phiMax_radm2=None,
                   verbose=False,
                   rm_verbose=False):
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
    # default connection (ie, local)
    with pymongo.MongoClient(host=host, connect=False) as dbclient:
        mydb = dbclient['spiceracs']  # Create/open database
        isl_col = mydb['islands']  # Create/open collection
        comp_col = mydb['components']  # Create/open collection
        beams_col = mydb['beams']  # Create/open collection

    # Basic querey
    myquery = {"Gaussian_ID": comp_id}
    doc = comp_col.find_one(myquery)

    iname = doc['Source_ID']
    cname = doc['Gaussian_ID']

    beams = beams_col.find_one({'Source_ID': iname})
    ifile = os.path.join(outdir, beams['beams'][field]['i_file'])
    outdir = os.path.dirname(ifile)

    header, dataI = do_RMsynth_3D.readFitsCube(ifile, rm_verbose)

    prefix = f'{outdir}/validation_{cname}'
    # Get source peak from Selavy
    ra = doc['RA']
    dec = doc['Dec']
    if len(dataI.shape) == 4:
        # drop Stokes axis
        wcs = WCS(header).dropaxis(2)
    else:
        wcs = WCS(header)

    x, y, z = np.array(wcs.all_world2pix(
        ra, dec, np.nanmean(freq), 0)).round().astype(int)

    mom = np.nansum(dataI, axis=0)

    plt.ion()
    plt.figure()
    plt.imshow(mom, origin='lower', cmap='cubehelix_r')
    plt.scatter(x, y, c='r', marker='x')
    plt.show()
    _ = input("Press [enter] to continue")  # wait for input from the user
    plt.close()    # close the figure to show the next one.

    data = np.nansum(dataI[:, y-1:y+1+1, x-1:x+1+1], axis=(1, 2))

    rmsi = estimate_noise_annulus(
        dataI.shape[2]//2,
        dataI.shape[1]//2,
        dataI
    )
    rmsi[rmsi == 0] = np.nan
    rmsi[np.isnan(rmsi)] = np.nanmedian(rmsi)
    noise = rmsi

    plt.ion()
    plt.figure()
    plt.step(freq/1e9, data)

    imod, qArr, uArr, dqArr, duArr, fitDict = \
        create_frac_spectra(freqArr=freq,
                            IArr=data,
                            QArr=data,
                            UArr=data,
                            dIArr=noise,
                            dQArr=noise,
                            dUArr=noise,
                            polyOrd=3,
                            verbose=True,
                            debug=False)
    plt.plot(freq/1e9, imod)
    plt.xlabel(r'$\nu$ [GHz]')
    plt.ylabel(r'$I$ [Jy/beam]')
    plt.tight_layout()
    plt.show()
    _ = input("Press [enter] to continue")  # wait for input from the user
    plt.close()    # close the figure to show the next one.

    data = data - imod
    data = data - np.nanmean(data)
    plt.ion()
    plt.figure()
    plt.step(freq/1e9, data)
    plt.xlabel(r'$\nu$ [GHz]')
    plt.ylabel(r'$I-\mathrm{model}(I)-\mathrm{mean}(\mathrm{model}(I))$')
    plt.tight_layout()
    plt.show()
    _ = input("Press [enter] to continue")  # wait for input from the user
    plt.close()    # close the figure to show the next one.

    datalist = [freq, data, data, dqArr, duArr]

    nSamples = nSamples
    phi_max = phiMax_radm2
    mDict, aDict = do_RMsynth_1D.run_rmsynth(
        datalist,
        phiMax_radm2=phi_max,
        nSamples=nSamples,
        verbose=True
    )
    plt.ion()
    plt.figure()
    plt.plot(aDict['phiArr_radm2'], abs(aDict['dirtyFDF']))
    plt.xlabel(r'$\phi$ [rad m$^{-2}$]')
    plt.ylabel(r'Dirty FDF (Stokes I)')
    plt.tight_layout()
    plt.show()
    _ = input("Press [enter] to continue")  # wait for input from the user
    plt.close()    # close the figure to show the next one.
    do_RMsynth_1D.saveOutput(mDict, aDict, prefix, verbose=verbose)


def main(field,
         outdir,
         host,
         client,
         dimension='1d',
         verbose=True,
         database=False,
         validate=False,
         limit=None,
         savePlots=False,
         weightType="variance",
         fitRMSF=True,
         phiMax_radm2=None,
         dPhi_radm2=None,
         nSamples=5,
         polyOrd=3,
         noStokesI=False,
         showPlots=False,
         not_RMSF=False,
         rm_verbose=False,
         debug=False,
         fit_function='log',
         tt0=None,
         tt1=None,
         ):

    outdir = os.path.abspath(outdir)
    outdir = os.path.join(outdir, 'cutouts')

    if savePlots:
        plotdir = os.path.join(outdir, 'plots')
        try_mkdir(plotdir)

    # default connection (ie, local)
    with pymongo.MongoClient(host=host, connect=False) as dbclient:
        mydb = dbclient['spiceracs']  # Create/open database
        isl_col = mydb['islands']  # Create/open collection
        comp_col = mydb['components']  # Create/open collection
        beams_col = mydb['beams']  # Create/open collection

    query = {
        '$and':  [
            {f'beams.{field}': {'$exists': True}},
            {f'beams.{field}.DR1': True}
        ]
    }

    beams = beams_col.find(query).sort('Source_ID')
    island_ids = sorted(beams_col.distinct('Source_ID', query))

    query = {'Source_ID': {'$in': island_ids}}
    components = comp_col.find(query).sort('Peak_flux', pymongo.DESCENDING)
    component_ids = [doc['Gaussian_ID'] for doc in components]

    n_comp = comp_col.count_documents(query)
    n_island = isl_col.count_documents(query)

    # Unset rmsynth in db
    if dimension == '1d':
        query = {
            '$and': [
                {
                    'Source_ID': {'$in': island_ids}
                },
                {
                    'rmsynth1d': True
                }
            ]
        }

        comp_col.update_many(query, {'$set': {'rmsynth1d': False}})

    elif dimension == '3d':
        query = {
            '$and': [
                {
                    'Source_ID': {'$in': island_ids}
                },
                {
                    'rmsynth3d': True
                }
            ]
        }

        isl_col.update(query, {'$set': {'rmsynth3d': False}})

    if limit is not None:
        count = limit
        n_comp = count
        n_island = count
        island_ids = island_ids[:count]
        component_ids = component_ids[:count]

    # Make frequency file
    freq, freqfile = getfreq(
        os.path.join(outdir, f"{beams[0]['beams'][f'{field}']['q_file']}"),
        outdir=outdir,
        filename='frequencies.txt',
        verbose=verbose
    )
    freq = np.array(freq)

    outputs = []

    if validate:
        if verbose:
            print(f'Running RMsynth on {n_comp} components')
        # We don't run this in parallel!
        for i, comp_id in enumerate(component_ids):
            output = rmsynthoncut_i(comp_id,
                                    outdir,
                                    freq,
                                    host,
                                    field,
                                    nSamples=nSamples,
                                    phiMax_radm2=phiMax_radm2,
                                    verbose=verbose,
                                    rm_verbose=rm_verbose)
            output.compute()

    elif dimension == '1d':
        if verbose:
            print(f'Running RMsynth on {n_comp} components')
        for i, comp_id in enumerate(component_ids):
            if i > n_comp+1:
                break
            else:
                output = rmsynthoncut1d(comp_id,
                                        outdir,
                                        freq,
                                        host,
                                        field,
                                        database=database,
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
                                        tt1=tt1
                                        )
                outputs.append(output)

    elif dimension == '3d':
        if verbose:
            print(f'Running RMsynth on {n_island} islands')

        for i, island_id in enumerate(island_ids):
            if i > n_island+1:
                break
            else:
                output = rmsynthoncut3d(island_id,
                                        outdir,
                                        freq,
                                        host,
                                        field,
                                        database=database,
                                        phiMax_radm2=phiMax_radm2,
                                        dPhi_radm2=dPhi_radm2,
                                        nSamples=nSamples,
                                        weightType=weightType,
                                        fitRMSF=fitRMSF,
                                        not_RMSF=not_RMSF,
                                        rm_verbose=rm_verbose
                                        )
                outputs.append(output)

    futures = client.persist(outputs)
    # dumb solution for https://github.com/dask/distributed/issues/4831
    time.sleep(5)
    tqdm_dask(futures, desc="Running RMsynth", disable=(not verbose))
    # progress(futures)

    if verbose:
        print('Done!')


def cli():
    """Command-line interface
    """
    import argparse
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
    from astropy.io.fits.verify import VerifyWarning
    warnings.simplefilter('ignore', category=VerifyWarning)
    warnings.simplefilter('ignore', category=RuntimeWarning)
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
    parser = argparse.ArgumentParser(description=descStr,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "field",
        metavar="field",
        type=str,
        help="RACS field to mosaic - e.g. 2132-50A."
    )
    parser.add_argument(
        'outdir',
        metavar='outdir',
        type=str,
        help='Directory containing cutouts (in subdir outdir/cutouts).')

    parser.add_argument(
        'host',
        metavar='host',
        type=str,
        help='Host of mongodb (probably $hostname -i).')

    parser.add_argument("--dimension", dest="dimension", default="1d",
                        help="How many dimensions for RMsynth [1d] or '3d'.")

    parser.add_argument("-v", dest="verbose", action="store_true",
                        help="verbose output [False].")

    parser.add_argument(
        "-m",
        dest="database",
        action="store_true",
        help="Add data to MongoDB [False]."
    )

    parser.add_argument(
        "--tt0",
        default=None,
        type=str,
        help="TT0 MFS image -- will be used for model of Stokes I -- also needs --tt1."
    )

    parser.add_argument(
        "--tt1",
        default=None,
        type=str,
        help="TT1 MFS image -- will be used for model of Stokes I -- also needs --tt0."
    )

    parser.add_argument("--validate", dest="validate", action="store_true",
                        help="Run on Stokes I [False].")

    parser.add_argument("--limit", dest="limit", default=None,
                        type=int, help="Limit number of sources [All].")

    # RM-tools args
    parser.add_argument("-sp", "--savePlots", action="store_true",
                        help="save the plots [False].")
    parser.add_argument("-w", dest="weightType", default="variance",
                        help="weighting [variance] (all 1s) or 'uniform'.")
    parser.add_argument("-f", dest="fit_function", type=str, default="log",
                        help="Stokes I fitting function: 'linear' or ['log'] polynomials.")
    parser.add_argument("-t", dest="fitRMSF", action="store_true",
                        help="Fit a Gaussian to the RMSF [False]")
    parser.add_argument("-l", dest="phiMax_radm2", type=float, default=None,
                        help="Absolute max Faraday depth sampled (overrides NSAMPLES) [Auto].")
    parser.add_argument("-d", dest="dPhi_radm2", type=float, default=None,
                        help="Width of Faraday depth channel [Auto].")
    parser.add_argument("-s", dest="nSamples", type=float, default=5,
                        help="Number of samples across the FWHM RMSF.")
    parser.add_argument("-o", dest="polyOrd", type=int, default=3,
                        help="polynomial order to fit to I spectrum [3].")
    parser.add_argument("-i", dest="noStokesI", action="store_true",
                        help="ignore the Stokes I spectrum [False].")
    parser.add_argument("-p", dest="showPlots", action="store_true",
                        help="show the plots [False].")
    parser.add_argument("-R", dest="not_RMSF", action="store_true",
                        help="Skip calculation of RMSF? [False]")
    parser.add_argument("-rmv", dest="rm_verbose", action="store_true",
                        help="Verbose RMsynth [False].")
    parser.add_argument("-D", dest="debug", action="store_true",
                        help="turn on debugging messages & plots [False].")

    args = parser.parse_args()

    if args.tt0 and not args.tt1:
        parser.error('the following arguments are required: tt1')
    elif args.tt1 and not args.tt0:
        parser.error('the following arguments are required: tt0')

    verbose = args.verbose

    cluster = LocalCluster(n_workers=12, processes=True,
                           threads_per_worker=1, local_directory='/dev/shm')
    client = Client(cluster)
    print(client)

    host = args.host
    if verbose:
        print('Testing MongoDB connection...')
    # default connection (ie, local)
    with pymongo.MongoClient(host=host, connect=False) as dbclient:
        try:
            dbclient.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")
        else:
            if verbose:
                print('MongoDB connection succesful!')

    main(field=args.field,
         outdir=args.outdir,
         host=host,
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
         tt1=args.tt1
         )


if __name__ == "__main__":
    cli()
