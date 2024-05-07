#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Make validation plots from a catalogue"""

import argparse
import logging
from pathlib import Path
from typing import NamedTuple as Struct

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import mad_std, sigma_clip
from astropy.table import Table
from astropy.wcs import WCS
from matplotlib.figure import Figure
from prefect import flow, task
from scipy import interpolate, stats

from arrakis.logger import UltimateHelpFormatter, logger
from arrakis.makecat import cat_parser
from arrakis.utils.pipeline import logo_str, upload_image_as_artifact_task
from arrakis.utils.typing import T

logger.setLevel(logging.INFO)


class GriddedMap(Struct):
    """Gridded catalogue data"""

    data: np.ndarray
    """Gridded data"""
    wcs: WCS
    """WCS of the gridded data"""


class BinnedMap(Struct):
    """Binned catalogue data"""

    data: np.ndarray
    """Binned data"""
    xc: np.ndarray
    """X bin centres"""
    yc: np.ndarray
    """Y bin centres"""
    wcs: WCS
    """WCS of the binned data"""


def make_gridded_map(
    tab: Table, column: str, npix: int = 512, map_size: u.Quantity = 8 * u.deg
) -> GriddedMap:
    logger.info(f"Making gridded map for {column}")
    coords = SkyCoord(ra=tab["ra"], dec=tab["dec"], unit="deg")
    coarse_shape = (npix, npix)
    coarse_header = fits.Header.fromstring(
        (
            f"NAXIS   =                    2\n"
            f"NAXIS1  =                  {coarse_shape[1]}\n"
            f"NAXIS2  =                  {coarse_shape[0]}\n"
            f"CTYPE1  = 'RA---SIN'\n"
            f"CRVAL1  = {coords.ra.deg.mean()}\n"
            f"CRPIX1  = {coarse_shape[1] / 2}\n"
            f"CDELT1  = {map_size.to(u.deg).value / coarse_shape[1]}\n"
            f"CUNIT1  = 'deg'\n"
            f"CTYPE2  = 'DEC--SIN'\n"
            f"CRVAL2  = {coords.dec.deg.mean()}\n"
            f"CRPIX2  = {coarse_shape[0] / 2}\n"
            f"CDELT2  = {map_size.to(u.deg).value / coarse_shape[0]}\n"
            f"CUNIT2  = 'deg'\n"
        ),
        sep="\n",
    )
    coarse_wcs = WCS(coarse_header)
    x, y = np.arange(coarse_shape[1]), np.arange(coarse_shape[0])
    X, Y = np.meshgrid(x, y)
    C = coarse_wcs.pixel_to_world(X, Y)
    R = C.ra.value
    D = C.dec.value
    sparse_points = np.stack(
        [coords.ra, coords.dec],
        -1,
    )  # shape (N, 2) in 2d

    data = interpolate.griddata(
        sparse_points,
        tab[column],
        (R, D),
        method="nearest",
    )  # default method is linear
    data_l = interpolate.griddata(
        sparse_points,
        tab[column],
        (R, D),
        method="linear",
    )  # default method is linear
    data[~np.isfinite(data_l)] = np.nan
    # Reverses the x-axis because RA
    data = data[:, ::-1]
    return GriddedMap(data, coarse_wcs)


def filter_then_median(arr: T) -> T:
    arr_clip = sigma_clip(
        arr, maxiters=None, sigma=3, cenfunc=np.nanmedian, stdfunc=mad_std
    )
    return np.nanmedian(arr_clip.compressed())


def make_binned_map(
    tab: Table,
    column: str,
    bins: int = 15,
    npix: int = 512,
    map_size: u.Quantity = 8 * u.deg,
) -> BinnedMap:
    logger.info(f"Making binned map for {column}")
    coords = SkyCoord(ra=tab["ra"], dec=tab["dec"], unit="deg")
    coarse_shape = (npix, npix)
    coarse_header = fits.Header.fromstring(
        (
            f"NAXIS   =                    2\n"
            f"NAXIS1  =                  {coarse_shape[1]}\n"
            f"NAXIS2  =                  {coarse_shape[0]}\n"
            f"CTYPE1  = 'RA---SIN'\n"
            f"CRVAL1  = {coords.ra.deg.mean()}\n"
            f"CRPIX1  = {coarse_shape[1] / 2}\n"
            f"CDELT1  = {map_size.to(u.deg).value / coarse_shape[1]}\n"
            f"CUNIT1  = 'deg'\n"
            f"CTYPE2  = 'DEC--SIN'\n"
            f"CRVAL2  = {coords.dec.deg.mean()}\n"
            f"CRPIX2  = {coarse_shape[0] / 2}\n"
            f"CDELT2  = {map_size.to(u.deg).value / coarse_shape[0]}\n"
            f"CUNIT2  = 'deg'\n"
        ),
        sep="\n",
    )
    coarse_wcs = WCS(coarse_header)
    x, y = coarse_wcs.world_to_pixel(coords)
    data, xe, ye, _ = stats.binned_statistic_2d(
        x, y, tab[column], statistic=filter_then_median, bins=bins
    )
    # Find bin centres
    xc = (xe[1:] + xe[:-1]) / 2
    yc = (ye[1:] + ye[:-1]) / 2
    XC, YC = np.meshgrid(xc, yc)
    # Reverses the x-axis because RA
    data = data[:, ::-1]
    return BinnedMap(data, XC, YC, coarse_wcs)


@task(name="rms and bkg plot")
def plot_rms_bkg(
    tab: Table,
    npix: int = 512,
    map_size: u.Quantity = 8 * u.deg,
) -> Figure:
    err_bkg_dict = {}
    for stokes in "IQU":
        err_bkg_dict[stokes] = {}
        for thing in ("err", "bkg"):
            err_bkg_dict[stokes][thing] = make_gridded_map(
                tab, f"stokes{stokes}_{thing}", npix=npix, map_size=map_size
            )
    mapping = {
        "I": ("I", "err"),
        "i": ("I", "bkg"),
        "Q": ("Q", "err"),
        "q": ("Q", "bkg"),
        "U": ("U", "err"),
        "u": ("U", "bkg"),
    }
    per_subplot_kw = {
        key: {"projection": err_bkg_dict[stokes][thing].wcs}
        for key, (stokes, thing) in mapping.items()
    }
    fig, ax_dict = plt.subplot_mosaic(
        """
        IQU
        iqu
        """,
        figsize=(24, 13),
        per_subplot_kw=per_subplot_kw,
        subplot_kw={
            "aspect": "equal",
        },
        sharex=True,
        sharey=True,
    )
    for key, ax in ax_dict.items():
        stokes, thing = mapping[key]
        data = err_bkg_dict[stokes][thing].data
        if thing == "err":
            im = ax.imshow(
                data * 1e6,
                origin="lower",
                cmap="YlOrRd",
                norm=plt.cm.colors.LogNorm(vmin=1e2, vmax=1e3),
            )
        else:
            im = ax.imshow(
                data * 1e6, origin="lower", cmap="coolwarm", vmin=-300, vmax=300
            )
        ax.set(xlabel="RA", ylabel="Dec")
        ax.grid()
        overlay = ax.get_coords_overlay("galactic")
        overlay.grid(color="tab:blue", ls="dashed", alpha=0.6)
        overlay[0].tick_params(colors="tab:blue")
        overlay[1].tick_params(colors="tab:blue")
        overlay[0].set_axislabel("$l$", color="tab:blue")
        overlay[1].set_axislabel("$b$", color="tab:blue")
        fig.colorbar(im, ax=ax, label="$\mu$Jy/beam", shrink=0.7, pad=0.15)
        ax.set_title(
            f"Stokes {stokes} {thing} - med: {np.nanmedian(data * 1e6)} $\mu$Jy/beam",
            pad=50,
        )

    return fig


@task(name="leakage plot")
def plot_leakage(
    tab: Table,
    snr_cut: float = 50,
    bins: int = 11,
    npix: int = 512,
    map_size: u.Quantity = 8 * u.deg,
) -> Figure:
    hi_i_tab = tab[tab["stokesI"] / tab["stokesI_err"] > snr_cut]
    hi_i_tab["stokesQ_frac"] = hi_i_tab["stokesQ"] / hi_i_tab["stokesI"]
    hi_i_tab["stokesU_frac"] = hi_i_tab["stokesU"] / hi_i_tab["stokesI"]
    leakage_dict = {}
    for stokes in "QU":
        leakage_dict[stokes] = make_binned_map(
            hi_i_tab, f"stokes{stokes}_frac", npix=npix, bins=bins, map_size=map_size
        )

    per_subplot_kw = {
        stokes: {"projection": val.wcs} for stokes, val in leakage_dict.items()
    }
    fig, ax_dict = plt.subplot_mosaic(
        """
        QU
        """,
        figsize=(16, 8),
        per_subplot_kw=per_subplot_kw,
        subplot_kw={
            "aspect": "equal",
        },
        sharex=True,
        sharey=True,
    )
    for stokes, ax in ax_dict.items():
        data = leakage_dict[stokes].data
        xc = leakage_dict[stokes].xc
        yc = leakage_dict[stokes].yc
        im = ax.pcolormesh(xc, yc, data, cmap="RdBu_r", vmin=-0.05, vmax=0.05)
        ax.set(xlabel="RA", ylabel="Dec")
        ax.grid()
        overlay = ax.get_coords_overlay("galactic")
        overlay.grid(color="tab:blue", ls="dashed", alpha=0.6)
        overlay[0].tick_params(colors="tab:blue")
        overlay[1].tick_params(colors="tab:blue")
        overlay[0].set_axislabel("$l$", color="tab:blue")
        overlay[1].set_axislabel("$b$", color="tab:blue")
        fig.colorbar(im, ax=ax, label="Fraction", shrink=0.7, pad=0.15)
        ax.set_title(
            f"Stokes {stokes}/I (binned) - absmed: {np.nanmedian(np.abs(data))*100}%",
            pad=50,
        )

    return fig


@flow(name="Validation")
def main(
    catalogue_path: Path,
    npix: int = 512,
    map_size: float = 8,
    snr_cut: float = 50,
    bins: int = 11,
):
    outdir = catalogue_path.parent
    tab = Table.read(catalogue_path)

    rms_bkg_fig = plot_rms_bkg(
        tab,
        npix=npix,
        map_size=map_size * u.deg,
    )
    rms_bkg_path = outdir / "validation_rms_bkg.png"
    rms_bkg_fig.savefig(rms_bkg_path)
    rms_bkg_uuid = upload_image_as_artifact_task(
        rms_bkg_path, description="Noise and background validation maps"
    )
    logger.info(f"Uploaded rms_bkg plot as {rms_bkg_uuid}")

    leakage_fig = plot_leakage(
        tab,
        snr_cut=snr_cut,
        bins=bins,
        npix=npix,
        map_size=map_size * u.deg,
    )
    leakage_path = outdir / "validation_leakage.png"
    leakage_fig.savefig(leakage_path)
    leakage_uuid = upload_image_as_artifact_task(
        leakage_path, description="Leakage validation maps"
    )
    logger.info(f"Uploaded leakage plot as {leakage_uuid}")

    logger.info("Validation plots complete")


def validation_parser(parent_parser: bool = False) -> argparse.ArgumentParser:
    descStr = f"""
    {logo_str}
    Arrakis:
    Validate RM catalogue.

    """
    val_parser = argparse.ArgumentParser(
        add_help=not parent_parser,
        description=descStr,
        formatter_class=UltimateHelpFormatter,
    )
    parser = val_parser.add_argument_group("validation options")
    parser.add_argument(
        "--npix",
        type=int,
        default=512,
        help="Number of pixels in the gridded maps",
    )
    parser.add_argument(
        "--map_size",
        type=float,
        default=8,
        help="Size of the maps in degrees",
    )
    return val_parser


def cli():
    catalogue_parser = cat_parser(parent_parser=True)
    val_parser = validation_parser(parent_parser=True)
    parser = argparse.ArgumentParser(
        parents=[val_parser, catalogue_parser],
        formatter_class=UltimateHelpFormatter,
        description=catalogue_parser.description,
    )
    args = parser.parse_args()

    main(
        catalogue_path=Path(args.outfile),
        npix=args.npix,
        map_size=args.map_size,
        snr_cut=args.leakage_snr,
        bins=args.leakage_bins,
    )


if __name__ == "__main__":
    cli()
