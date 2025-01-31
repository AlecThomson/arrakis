#!/usr/bin/env python
"""Make validation plots from a catalogue"""

from __future__ import annotations

import argparse
import logging
from importlib import resources
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
    try:
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
                f"Stokes {stokes} {thing} - med: {np.nanmedian(data * 1e6):0.1f}$\pm${np.nanstd(data * 1e6):0.1f} $\mu$Jy/beam\n - min: {np.nanmin(data * 1e6):0.1f} $\mu$Jy/beam",
                pad=50,
            )

        return fig
    except Exception as e:
        logger.error(f"Error in making rms and bkg plot: {e}")
        return plt.figure()


@task(name="leakage plot")
def plot_leakage(
    tab: Table,
    snr_cut: float = 50,
    bins: int = 11,
    npix: int = 512,
    map_size: u.Quantity = 8 * u.deg,
) -> Figure:
    try:
        hi_i_tab = tab[tab["stokesI"] / tab["stokesI_err"] > snr_cut]
        hi_i_tab["stokesQ_frac"] = hi_i_tab["stokesQ"] / hi_i_tab["stokesI"]
        hi_i_tab["stokesU_frac"] = hi_i_tab["stokesU"] / hi_i_tab["stokesI"]
        leakage_dict = {}
        for stokes in "QU":
            try:
                leakage_dict[stokes] = make_binned_map(
                    hi_i_tab,
                    f"stokes{stokes}_frac",
                    npix=npix,
                    bins=bins,
                    map_size=map_size,
                )
            except Exception as e:
                logger.error(f"Error in making binned map for {stokes}: {e}")
                leakage_dict[stokes] = BinnedMap(
                    np.full((bins, bins), np.nan),
                    np.linspace(-4, 4, bins),
                    np.linspace(-4, 4, bins),
                    None,
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
            if leakage_dict[stokes].wcs is None:
                continue
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
                f"Stokes {stokes}/I (binned) - absmed: {np.nanmedian(np.abs(data)) * 100:0.1f}$\pm${np.nanstd(np.abs(data)) * 100:0.1f}%",
                pad=50,
            )

        return fig
    except Exception as e:
        logger.error(f"Error in making leakage plot: {e}")
        return plt.figure()


def cross_match(
    my_tab: Table, other_tab: Table, radius: u.Quantity = 1 * u.arcsec
) -> Table:
    my_coords = SkyCoord(ra=my_tab["ra"], dec=my_tab["dec"], unit="deg")
    other_coords = SkyCoord(ra=other_tab["ra"], dec=other_tab["dec"], unit="deg")
    idx, d2d, _ = my_coords.match_to_catalog_sky(other_coords)
    sep_constraint = d2d < radius
    other_match = other_tab[idx[sep_constraint]]
    my_match = my_tab[sep_constraint]
    return my_match, other_match


@task(name="rm plot")
def plot_rm(
    tab: Table,
    npix: int = 512,
    map_size: u.Quantity = 8 * u.deg,
) -> Figure:
    try:
        good_idx = (
            (~tab["snr_flag"])
            & (~tab["leakage_flag"])
            & (~tab["channel_flag"])
            & (~tab["stokesI_fit_flag"])
            & (~tab["local_rm_flag"])
        )

        good_tab = tab[good_idx]

        nvss_path = resources.files("arrakis.data") / "Taylor2009.fits.zip"
        nvss_tab = Table.read(nvss_path, format="fits")
        spass_path = resources.files("arrakis.data") / "Schnitzeler2019.fits.zip"
        spass_tab = Table.read(spass_path, format="fits")

        rm_gridded = make_gridded_map(good_tab, "rm", npix=npix, map_size=map_size)

        fig, ax_dict = plt.subplot_mosaic(
            """
            MNS
            """,
            figsize=(16, 4),
            per_subplot_kw={
                "M": {"projection": rm_gridded.wcs},
            },
        )
        for label, other_cat, ax in zip(
            ("NVSS", "SPASS"),
            (nvss_tab, spass_tab),
            (ax_dict["N"], ax_dict["S"]),
        ):
            ax.set_title(label)
            racs_match, other_match = cross_match(
                good_tab, other_cat, radius=60 * u.arcsec
            )
            if len(racs_match) == 0:
                # Hide axes
                ax.axis("off")
                ax.text(
                    x=0.5,
                    y=0.5,
                    s=f"No {label} matches",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                )
                continue
            _ = ax.errorbar(
                racs_match["rm"],
                other_match["rm"],
                xerr=racs_match["rm_err"] * 5,
                yerr=other_match["rm_err"] * 5,
                fmt="o",
                label="$\pm 5 \sigma$",
            )
            abs_max_val = np.nanmax(
                np.abs(np.concatenate([racs_match["rm"], other_match["rm"]]))
            )
            ax.plot(
                [-abs_max_val, abs_max_val],
                [-abs_max_val, abs_max_val],
                color="k",
                ls="--",
            )
            ax.legend()

            ax.set(
                xlabel=f"RACS  / {u.rad / u.m**2:latex_inline}",
                ylabel=f"{label} / {u.rad / u.m**2:latex_inline}",
                aspect="equal",
            )
        _ = ax_dict["M"].imshow(
            rm_gridded.data * np.nan,
            origin="lower",
            cmap="coolwarm",
            vmin=-100,
            vmax=100,
        )

        def rm_scaler(rm: np.ndarray) -> np.ndarray:
            # Scale the RM value by 2
            return 2 * np.abs(rm)

        pos_idx = good_tab["rm"] > 0
        neg_idx = good_tab["rm"] < 0

        _ = ax_dict["M"].scatter(
            good_tab["ra"][pos_idx],
            good_tab["dec"][pos_idx],
            edgecolor="tab:red",
            s=rm_scaler(good_tab["rm"][pos_idx]),
            facecolor="none",
            transform=ax_dict["M"].get_transform("world"),
            linewidths=1,
        )

        _ = ax_dict["M"].scatter(
            good_tab["ra"][neg_idx],
            good_tab["dec"][neg_idx],
            edgecolor="tab:blue",
            s=rm_scaler(good_tab["rm"][neg_idx]),
            facecolor="none",
            transform=ax_dict["M"].get_transform("world"),
            linewidths=1,
        )

        for rm in [-100, -10, +10, 100]:
            _ = ax_dict["M"].scatter(
                np.nan,
                np.nan,
                edgecolor="tab:blue" if rm < 0 else "tab:red",
                s=rm_scaler(rm),
                facecolor="none",
                transform=ax_dict["M"].get_transform("world"),
                label=rf"{rm}",
                linewidths=1,
            )

        ax_dict["M"].legend()

        ax_dict["M"].set(xlabel="RA", ylabel="Dec")
        ax_dict["M"].grid()
        overlay = ax_dict["M"].get_coords_overlay("galactic")
        overlay.grid(color="tab:blue", ls="dashed", alpha=0.6)
        overlay[0].tick_params(colors="tab:blue")
        overlay[1].tick_params(colors="tab:blue")
        overlay[0].set_axislabel("$l$", color="tab:blue")
        overlay[1].set_axislabel("$b$", color="tab:blue")
        ax_dict["M"].set_title("RM bubble", pad=50)
        fig.suptitle(
            "rotation measure",
            y=1.1,
        )
        return fig
    except Exception as e:
        logger.error(f"Error in making RM plot: {e}")
        return plt.figure()


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

    try:
        rms_bkg_fig = plot_rms_bkg(
            tab,
            npix=npix,
            map_size=map_size * u.deg,
        )
        rms_bkg_path = outdir / "validation_rms_bkg.png"
        rms_bkg_fig.savefig(rms_bkg_path, bbox_inches="tight", dpi=72)
        rms_bkg_uuid = upload_image_as_artifact_task(
            rms_bkg_path, description="Noise and background validation maps"
        )
        logger.info(f"Uploaded rms_bkg plot as {rms_bkg_uuid}")
    except Exception as e:
        logger.error(f"Error in rms_bkg plot: {e}")

    try:
        leakage_fig = plot_leakage(
            tab,
            snr_cut=snr_cut,
            bins=bins,
            npix=npix,
            map_size=map_size * u.deg,
        )
        leakage_path = outdir / "validation_leakage.png"
        leakage_fig.savefig(leakage_path, bbox_inches="tight", dpi=72)
        leakage_uuid = upload_image_as_artifact_task(
            leakage_path, description="Leakage validation maps"
        )
        logger.info(f"Uploaded leakage plot as {leakage_uuid}")
    except Exception as e:
        logger.error(f"Error in leakage plot: {e}")

    try:
        rm_fig = plot_rm(
            tab,
            npix=npix,
            map_size=map_size * u.deg,
        )
        rm_path = outdir / "validation_rm.png"
        rm_fig.savefig(rm_path, bbox_inches="tight", dpi=72)
        rm_uuid = upload_image_as_artifact_task(
            rm_path, description="Rotation measure validation maps"
        )
        logger.info(f"Uploaded rm plot as {rm_uuid}")
    except Exception as e:
        logger.error(f"Error in rm plot: {e}")

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
