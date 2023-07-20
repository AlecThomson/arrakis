#!/usr/bin/env python3
import json

import astropy
import astropy.units as units
import matplotlib.pyplot as plt
import numpy as np
import pymongo
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from tqdm.auto import tqdm, trange

from arrakis.logger import logger, logging
from arrakis.utils.database import get_db
from arrakis.utils.fitsutils import getdata
from arrakis.utils.json import MyEncoder


def makesurf(start, stop, field, datadir, save_plots=True, data=None):
    # myquery = {'rmsynth1d': True}
    query = {"$and": [{f"beams.{field}": {"$exists": True}}]}

    beams = list(beams_col.find(query).sort("Source_ID"))
    island_ids = sorted(beams_col.distinct("Source_ID", query))

    query = {"Source_ID": {"$in": island_ids}}
    # myquery = {'rmsynth1d': True}
    components = list(comp_col.find(query).sort("Source_ID"))
    ras, decs, freqs, stokeis, stokeqs, stokeus = [], [], [], [], [], []
    specs = []
    for i, comp in enumerate(tqdm(components)):
        iname = comp["Source_ID"]
        cname = comp["Gaussian_ID"]
        spectra = f"{datadir}/cutouts/{iname}/{cname}.dat"
        if data is None:
            try:
                freq, iarr, qarr, uarr, rmsi, rmsq, rmsu = np.loadtxt(spectra).T
                specs.append([freq, iarr, qarr, uarr, rmsi, rmsq, rmsu])
            except Exception as e:
                logger.warning(f"Could not find '{spectra}': {e}")
                continue
        else:
            try:
                freq, iarr, qarr, uarr, rmsi, rmsq, rmsu = data[i]
            except IndexError:
                continue
        ras.append(comp["RA"])
        decs.append(comp["Dec"])
        freqs.append(np.nanmean(freq[start:stop]))
        stokeis.append(np.nansum(iarr[start:stop]))
        stokeqs.append(np.nansum(qarr[start:stop]))
        stokeus.append(np.nansum(uarr[start:stop]))

    ras = np.array(ras)
    decs = np.array(decs)
    stokeis = np.array(stokeis)
    stokeqs = np.array(stokeqs)
    stokeus = np.array(stokeus)
    freqs = np.nanmean(np.array(freqs))
    logger.debug("freq is ", freqs)
    coords = SkyCoord(ras * units.deg, decs * units.deg)
    wcs = WCS(
        f"/group/askap/athomson/projects/RACS/CI0_mosaic_1.0/RACS_test4_1.05_{field}.fits"
    )
    x, y = wcs.celestial.world_to_pixel(coords)
    # Parse out data
    x_raw = x  # ra
    y_raw = y  # dec
    q_raw = stokeqs / stokeis  # Q/I
    u_raw = stokeus / stokeis  # U/I

    # Kill nans
    good_idxs = (~np.isnan(q_raw)) & (~np.isnan(u_raw))

    # Use good data only
    x = x_raw[good_idxs]
    y = y_raw[good_idxs]
    q = q_raw[good_idxs]
    u = u_raw[good_idxs]

    # x = x/60**2 #get emil's pixel positions in degs
    # y = y/60**2 #get emil's pixel positions in degs
    p = np.sqrt(q.astype(float) ** 2 + u.astype(float) ** 2)

    # Imports
    from scipy.spatial import distance_matrix

    # Settings
    pixelscales = astropy.wcs.utils.proj_plane_pixel_scales(wcs)
    dperpix = pixelscales[1]
    d = (1 / 60) * 15 / dperpix  # Radius for circular estimator sliding window aperture
    trim_mean_frac = (
        0.2  # frac of data points to chop from each end of frac Stokes value dist
    )
    # Minimum number of sources required in the sliding aperture to return a non-nan estimate of the local leakage
    min_data_points_in_aperture = 5
    grid_point_sep_deg = d / 4

    # Define functions
    def trim_mean(x):
        from scipy import stats

        return stats.trim_mean(x, trim_mean_frac)

    # Positions of measured leakages
    pos_measurements = np.array(list(zip(x, y)))

    # Positions of grid points to derive leakage estimates at
    xnew = np.arange(np.min(x), np.max(x) + grid_point_sep_deg, grid_point_sep_deg)
    ynew = np.arange(np.min(y), np.max(y) + grid_point_sep_deg, grid_point_sep_deg)
    logger.debug(len(xnew), len(ynew))
    xxnew, yynew = np.meshgrid(xnew, ynew)
    pos_estimator_grid = np.array([[a, b] for a in xnew for b in ynew])

    # Calculate pair-wise distances between the two sets of coordinate pairs
    logger.info("\nDeriving pair-wise distance matrix...")
    pair_dist = distance_matrix(pos_estimator_grid, pos_measurements)
    logger.info("Done.\n")

    # Collect leakage values nearby each grid point
    q_estimates = []
    u_estimates = []
    p_estimates = []
    num_points_in_aperture_list = []  # Init collectors

    logger.info("\nDeriving robust leakage estimates for interpolation grid...")
    for row_idx, row in enumerate(tqdm(pair_dist)):
        # Guide to where we're at
        # if row_idx%100==0:
        # 	logger.info('Processing row %d of %d'%(row_idx,len(pair_dist)))

        # idxs of poitns within d degs
        idxs_of_points_in_aperture = np.argwhere(row < d)
        # collect data points for sources in aperture
        q_of_points_in_aperture = q[idxs_of_points_in_aperture]
        u_of_points_in_aperture = u[idxs_of_points_in_aperture]
        p_of_points_in_aperture = p[idxs_of_points_in_aperture]
        # robust estimator of central value of dist
        if len(q_of_points_in_aperture) >= min_data_points_in_aperture:
            est_q_leak_of_points_in_aperture = trim_mean(q_of_points_in_aperture)
            q_estimates.append(est_q_leak_of_points_in_aperture)
            num_points_in_aperture_list.append(len(q_of_points_in_aperture))
        else:
            q_estimates.append(np.nan)
        if len(u_of_points_in_aperture) >= min_data_points_in_aperture:
            est_u_leak_of_points_in_aperture = trim_mean(u_of_points_in_aperture)
            u_estimates.append(est_u_leak_of_points_in_aperture)
        else:
            u_estimates.append(np.nan)
        if len(p_of_points_in_aperture) >= min_data_points_in_aperture:
            est_p_leak_of_points_in_aperture = trim_mean(p_of_points_in_aperture)
            p_estimates.append(est_p_leak_of_points_in_aperture)
        else:
            p_estimates.append(np.nan)

    q_estimates_arr = np.array(q_estimates)
    u_estimates_arr = np.array(u_estimates)
    p_estimates_arr = np.array(p_estimates)

    logger.info(
        "\nThe mean number of points in each aperture of %.2f degs was %d\n"
        % (d, np.nanmean(num_points_in_aperture_list))
    )

    # plot results
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(
        111,
    )
    # q_leakage_map = np.rot90(q_estimates_arr.reshape((len(xnew),len(ynew))).astype(float),k=3)
    q_leakage_map = np.rot90(
        q_estimates_arr.reshape((len(xnew), len(ynew))).astype(float), k=3
    )
    im = ax.imshow(
        q_leakage_map, origin="lower", vmin=-0.05, vmax=0.05, cmap="coolwarm"
    )
    fig.colorbar(im, label="Q/I")
    ax.set_aspect("equal", "box")
    ax.invert_xaxis()
    if save_plots:
        plt.savefig(f"q_leakage_{field}.png")
    # plt.xlim(-10,60)
    # plt.ylim(-10,40)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(
        111,
    )
    u_leakage_map = np.rot90(
        u_estimates_arr.reshape((len(xnew), len(ynew))).astype(float), k=3
    )
    im = ax.imshow(
        u_leakage_map, origin="lower", vmin=-0.05, vmax=0.05, cmap="coolwarm"
    )
    fig.colorbar(im, label="U/I")
    ax.set_aspect("equal", "box")
    ax.invert_xaxis()
    if save_plots:
        plt.savefig(f"u_leakage_{field}.png")
    # plt.xlim(-10,60)
    # plt.ylim(-10,40)
    return freqs, q_leakage_map, u_leakage_map, specs, wcs


def main(field, datadir, username="admin", password=None):
    global beams_col
    global island_col
    global comp_col

    beams_col, island_col, comp_col = get_db(
        "146.118.68.63", username=username, password=password
    )

    start = 0
    stop = -1
    freqs, q_leakage_map, u_leakage_map, data, wcs = makesurf(
        start, stop, field, datadir, save_plots=True
    )

    start = 0
    f_big, q_big, u_big = [], [], []
    for i in trange(6):
        stop = start + (288 // 6) - 1
        f, q, u, _, _ = makesurf(
            start, stop, field, datadir, save_plots=False, data=data
        )
        f_big.append(f)
        q_big.append(q)
        u_big.append(u)
        start += (288 // 6) - 1
    f_big = np.array(f_big)
    q_big = np.array(q_big)
    u_big = np.array(u_big)

    fig, ax = plt.subplots(
        2,
        6,
        figsize=(18, 6),
    )
    lim = 0.1
    for i, (f, q, u) in enumerate(zip(f_big, q_big, u_big)):
        ax[0, i].imshow(q, origin="lower", vmin=-lim, vmax=lim, cmap=plt.cm.coolwarm)
        ax[0, i].axis("off")
        ax[0, i].invert_xaxis()
        ax[0, i].set_title(f"{f/1e6:0.1f}MHz")
        ax[1, i].imshow(
            u,
            origin="lower",
            cmap=plt.cm.coolwarm,
            vmin=-lim,
            vmax=lim,
        )
        ax[1, i].axis("off")
        ax[1, i].invert_xaxis()

    ax[0, 0].text(-20, 50, "Q")
    ax[1, 0].text(-20, 50, "U")
    # plt.subplots_adjust(hspace=0)
    plt.savefig(f"{field}_leakages.png")


def cli():
    import argparse
    import getpass

    # Help string to be shown using the -h option
    descStr = f"""
    Make leakage plots for a field.

    """
    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field to mosaic - e.g. 2132-50A."
    )

    parser.add_argument(
        "datadir",
        metavar="datadir",
        type=str,
        help="Directory containing cutouts (in subdir outdir/cutouts)..",
    )

    args = parser.parse_args()
    password = getpass.getpass()
    main(field=args.field, datadir=args.datadir, username="admin", password=password)


if __name__ == "__main__":
    cli()
