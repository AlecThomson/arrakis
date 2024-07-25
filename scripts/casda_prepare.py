#!/usr/bin/env python3
"""Prepare files for CASDA upload"""

from __future__ import annotations

import argparse
import hashlib
import logging
import os
import pickle
import subprocess as sp
import tarfile
import time
from glob import glob

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polspectra
from arrakis.logger import TqdmToLogger, logger
from arrakis.makecat import write_votable
from arrakis.utils.io import try_mkdir, try_symlink
from arrakis.utils.pipeline import chunk_dask
from astropy.io import fits
from astropy.table import Column, Table
from astropy.units.core import get_current_unit_registry
from astropy.visualization import ImageNormalize, MinMaxInterval, SqrtStretch
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from dask import delayed
from prefect_dask import get_dask_client
from radio_beam import Beam
from spectral_cube.cube_utils import convert_bunit
from tqdm.auto import tqdm

TQDM_OUT = TqdmToLogger(logger, level=logging.INFO)


def make_thumbnail(cube_f: str, cube_dir: str):
    cube = fits.getdata(cube_f)
    _ = cube[:, 0, :, :]
    qcube = cube[:, 1, :, :]
    ucube = cube[:, 2, :, :]
    pcube = np.hypot(qcube, ucube)

    head = fits.getheader(cube_f)
    wcs = WCS(head)
    med = np.nanmedian(pcube, axis=0)
    med_wcs = wcs.celestial
    pix_scales = proj_plane_pixel_scales(med_wcs) * u.deg
    beam = Beam.from_fits_header(head)
    ellipse = beam.ellipse_to_plot(
        xcen=10,
        ycen=10,
        pixscale=pix_scales[0],  # Assume same pixel scale in x and y
    )
    fig = plt.figure(facecolor="w")
    cmap = "gray" if ".weights." in cube_f else "Purples"
    norm = ImageNormalize(med, interval=MinMaxInterval(), stretch=SqrtStretch())
    ax = fig.add_subplot(111, projection=med_wcs)
    im = ax.imshow(med, origin="lower", cmap=cmap, norm=norm)
    ax.add_patch(ellipse)
    ellipse.set_clip_box(ax.bbox)
    ellipse.set_facecolor("none")
    ellipse.set_edgecolor("k")
    ellipse.set_hatch("//")
    ax.set_xlabel("RA")
    ax.set_ylabel("Dec")
    fig.colorbar(im, ax=ax, label=f"{convert_bunit(head['BUNIT']):latex_inline}")
    outf = os.path.join(cube_dir, os.path.basename(cube_f).replace(".fits", ".png"))
    logger.info(f"Saving thumbnail to {outf}")
    fig.savefig(outf, dpi=72)
    plt.close(fig)


def find_spectra(data_dir: str = ".") -> list:
    """Find spectra in from cutouts directory

    Args:
        data_dir (str, optional): Directory containing cutouts directory. Defaults to ".".

    Returns:
        list: List of spectra in ascii format
    """
    cut_dir = os.path.join(data_dir, "cutouts")
    logger.info(f"Globbing for spectra in {cut_dir}")
    spectra = glob(os.path.join(os.path.join(cut_dir, "*"), "*[0-9].dat"))
    logger.info(f"Found {len(spectra)} spectra (in frequency space)")
    return spectra


@delayed()
def convert_spectra(
    number: int,
    spectrum: str,
    fname_polcat_hash: str,
    spec_dir: str = ".",
) -> str:
    """Convert a ascii spectrum to FITS

    Args:
        spectrum (str): Name of ASCII spectrum file
        spec_dir (str, optional): Directory to save FITS spectrum. Defaults to '.'.
    """
    # Re-register astropy units
    for unit in (u.deg, u.hour, u.hourangle, u.Jy, u.arcsec, u.arcmin, u.beam):
        get_current_unit_registry().add_enabled_units([unit])
    with open(fname_polcat_hash, "rb") as f:
        cat_row = Table(pickle.load(f)[number])
    rmsf_unit = u.def_unit("RMSF")
    unit_flux = u.Jy / u.beam
    unit_fdf = unit_flux / rmsf_unit
    radms = u.radian / u.m**2
    # First deal with the frequency data
    full_data = pd.read_csv(
        spectrum,
        names=("freq", "I", "Q", "U", "dI", "dQ", "dU"),
        delim_whitespace=True,
    )

    freq = full_data["freq"].values
    u.deg = u.core._recreate_irreducible_unit(u.Unit, ["deg", "degree"], True)
    spectrum_table = polspectra.from_arrays(
        long_array=cat_row["ra"].value,
        lat_array=cat_row["dec"].value,
        freq_array=freq,
        stokesI=(full_data.I.values)[np.newaxis, :],
        stokesI_error=(full_data.dI.values)[np.newaxis, :],
        stokesQ=(full_data.Q.values)[np.newaxis, :],
        stokesQ_error=(full_data.dQ.values)[np.newaxis, :],
        stokesU=(full_data.U.values)[np.newaxis, :],
        stokesU_error=(full_data.dU.values)[np.newaxis, :],
        source_number_array=[number],
        cat_id=cat_row["cat_id"],
        beam_maj=cat_row["beam_maj"],
        beam_min=cat_row["beam_min"],
        beam_pa=cat_row["beam_pa"],
        coordinate_system="icrs",
        channel_width=np.diff(freq)[0],
        telescope=cat_row["telescope"],
        epoch=cat_row["epoch"],
        integration_time=cat_row["int_time"],
        leakage=cat_row["leakage"],
        flux_type=cat_row["flux_type"],
    )
    # Fix units
    for col in (
        "stokesI",
        "stokesI_error",
        "stokesQ",
        "stokesQ_error",
        "stokesU",
        "stokesU_error",
    ):
        spectrum_table[col].unit = unit_flux

    pol_df_cols = {
        "faraday_depth": {
            "unit": radms,
            "description": "Faraday depth",
        },
        "faraday_depth_long": {
            "unit": radms,
            "description": "Faraday depth (long)",
        },
        "FDF_Q_dirty": {
            "unit": unit_fdf,
            "description": "Dirty Stokes Q per Faraday depth",
        },
        "FDF_U_dirty": {
            "unit": unit_fdf,
            "description": "Dirty Stokes U per Faraday depth",
        },
        "FDF_Q_clean": {
            "unit": unit_fdf,
            "description": "Clean Stokes Q per Faraday depth",
        },
        "FDF_U_clean": {
            "unit": unit_fdf,
            "description": "Clean Stokes U per Faraday depth",
        },
        "FDF_Q_model": {
            "unit": unit_fdf,
            "description": "Model Stokes Q per Faraday depth",
        },
        "FDF_U_model": {
            "unit": unit_fdf,
            "description": "Model Stokes U per Faraday depth",
        },
        "RMSF_Q": {
            "unit": unit_fdf,
            "description": "Stokes Q RMSF per Faraday depth",
        },
        "RMSF_U": {
            "unit": unit_fdf,
            "description": "Stokes U RMSF per Faraday depth",
        },
    }

    # Now deal with the RM data
    suffixes = (
        "_FDFclean.dat",
        "_FDFdirty.dat",
        "_FDFmodel.dat",
        "_RMSF.dat",
    )
    rmtables = {}  # type: Dict[str, Table]
    for suffix in suffixes:
        name = suffix.replace(".dat", "").replace("_", "").replace("FDF", "")
        rm_file = spectrum.replace(".dat", suffix)
        full_rm_data = pd.read_csv(
            rm_file, names=("phi", "Q", "U"), delim_whitespace=True
        )
        rmtables[name] = full_rm_data

    for col, desc in tqdm(
        pol_df_cols.items(), desc="Adding spectrum columns", file=TQDM_OUT
    ):
        keys = col.split("_")
        if keys[0] == "faraday":
            if keys[-1] == "long":
                data = rmtables["RMSF"]["phi"].values * radms
            else:
                data = rmtables["clean"]["phi"].values * radms
        elif keys[0] == "RMSF":
            data = rmtables["RMSF"][keys[1]].values * u.Jy / u.beam / rmsf_unit
        elif keys[0] == "FDF":
            data = rmtables[keys[2]][keys[1]].values * u.Jy / u.beam / rmsf_unit
        else:
            msg = f"Unknown column {col}"
            raise ValueError(msg)

        new_col = Column(
            name=col,
            unit=data.unit,
            dtype="object",
            shape=(),
            length=spectrum_table.Nrows,
            description=desc["description"],
        )
        new_col[:] = [data.value]
        spectrum_table.table.add_column(new_col)

    outf = os.path.join(
        spec_dir, os.path.basename(spectrum).replace(".dat", "_polspec.fits")
    )
    logger.info(f"Saving to {outf}")
    spectrum_table.write_FITS(outf, overwrite=True)
    # Add object to header
    # Hard code the pixel size for now
    pix_size = 2.5 * u.arcsec
    with fits.open(outf, mode="update") as hdul:
        hdul[0].header["OBJECT"] = (cat_row["cat_id"][0], "Gaussian ID")
        hdul[0].header["NAXIS"] = 2
        hdul[0].header["NAXIS1"] = 1
        hdul[0].header["NAXIS2"] = 1
        hdul[0].header["CTYPE1"] = "RA---SIN"
        hdul[0].header["CTYPE2"] = "DEC--SIN"
        hdul[0].header["CRVAL1"] = cat_row["ra"][0]
        hdul[0].header["CRVAL2"] = cat_row["dec"][0]
        hdul[0].header["CDELT1"] = pix_size.to(u.deg).value
        hdul[0].header["CDELT2"] = pix_size.to(u.deg).value
        hdul[0].header["CRPIX1"] = 1
        hdul[0].header["CRPIX2"] = 1
        hdul[0].header["CUNTI1"] = "deg"
        hdul[0].header["CUNTI2"] = "deg"
        hdul[0].header["comment"] = (
            "Dummy image to indicate the pixel size and position"
        )
        # Add dummy data to make it a valid FITS file
        hdul[0].data = np.zeros((1, 1))
        hdul.flush()

    return outf


@delayed()
def update_cube(cube: str, cube_dir: str) -> None:
    """Update cube headers and symlink to CASDA area

    Args:
        cube (str): Cubelet path
        cube_dir (str): CASDA cublet directory
    """
    # Re-register astropy units
    for unit in (u.deg, u.hour, u.hourangle, u.Jy, u.arcsec, u.arcmin, u.beam):
        get_current_unit_registry().add_enabled_units([unit])
    stokes = ("i", "q", "u")
    imtypes = ("image.restored", "weights")
    idata = fits.getdata(
        cube,
        memmap=True,
    )
    cube = os.path.abspath(cube)

    for imtype in imtypes:
        data = np.zeros((idata.shape[0], len(stokes), idata.shape[2], idata.shape[3]))
        for i, stoke in enumerate(stokes):
            fname = cube.replace("image.restored.i", f"{imtype}.{stoke}")
            if stoke != "i":
                fname = fname.replace(
                    ".linmos.edge.linmos.fits", ".linmos.ion.edge.linmos.fits"
                )
                if not os.path.exists(fname) and imtype == "weights":
                    fname = fname.replace(
                        ".linmos.ion.edge.linmos.fits", ".linmos.edge.linmos.fits"
                    )
            if not os.path.exists(fname):
                msg = f"Could not find {fname}"
                raise FileNotFoundError(msg)
            data[:, i, :, :] = fits.getdata(
                fname,
                memmap=True,
            )[:, 0, :, :]

            # Get header from Stokes Q
            if stoke == "q":
                header = fits.getheader(fname)
        # Get source ID and update header
        source_id = os.path.basename(os.path.dirname(cube))
        header["OBJECT"] = (source_id, "Source ID")
        header["CRVAL3"] = 1.0
        if header["BUNIT"] == "beam-1 Jy":
            header["BUNIT"] = (u.Jy / u.beam).to_string()

        outf = os.path.join(
            cube_dir,
            os.path.basename(cube).replace(
                "image.restored.i.", f"{imtype}.{''.join(stokes)}."
            ),
        ).replace("RACS_test4_1.05_", "RACS_")
        logger.info(f"Writing {outf} cubelet")
        fits.writeto(outf, data, header, overwrite=True)

        # Move cube to cubelets directory
        make_thumbnail(outf, cube_dir)


def find_cubes(data_dir: str = ".") -> list:
    """Find Stokes I image cubelets in a directory

    Args:
        data_dir (str, optional): Data containg cutouts directory. Defaults to ".".

    Returns:
        list: List of cubelets
    """
    cut_dir = os.path.join(data_dir, "cutouts")
    logger.info(f"Globbing for cubes in {cut_dir}")
    cubes = glob(
        os.path.join(os.path.join(cut_dir, "*"), "*.image.restored.i.*.linmos.fits")
    )
    logger.info(f"Found {len(cubes)} Stokes I image cubes")
    return cubes


# @delayed(nout=2)
def init_polspec(
    casda_dir: str,
    spectrum_table_0: str,
    outdir: str = "",
):
    if outdir == "":
        outdir = casda_dir
    polspec_0 = polspectra.from_FITS(spectrum_table_0)
    out_fits = os.path.join(os.path.abspath(outdir), "spice_racs_dr1_polspec.fits")
    logger.info(f"Saving to {out_fits}")
    polspec_0.write_FITS(out_fits, overwrite=True)

    out_hdf = os.path.join(os.path.abspath(outdir), "spice_racs_dr1_polspec.hdf5")
    # logger.info(f"Saving to {out_hdf}")
    # polspec_0.write_HDF5(out_hdf, overwrite=True, compress=True)

    return out_fits, out_hdf


# Reworked from polspectra.py
def write_polspec(table: Table, filename: str, overwrite: bool = False):
    """Write the polspectra to a FITS file.
    Parameters:
        table: astropy.table.Table
        filename : str
            Name and relative path of the file to save to.
        overwrite : bool [False]
    Overwrite the file if it already exists?
    """
    # This is going to be complicated, because the automatic write algorithm
    # doesn't like variable length arrays. pyfits can support it, it just
    # needs a little TLC to get it into the correct format.
    # Converting to FITSrec format loses the column description...
    # Is that important?
    fits_columns = []
    col_descriptions = []

    # per column, convert to fits column:
    for col in table.colnames:
        tabcol = table[col]
        if tabcol.dtype != np.dtype("object"):  # Normal columns
            col_format = fits.column._convert_record2fits(tabcol.dtype)
        else:  # Channelized columns
            subtype = np.result_type(
                np.array(tabcol[0])
            )  # get the type of each element in 2D array
            col_format = "Q" + fits.column._convert_record2fits(subtype) + "()"
        unit = tabcol.unit.to_string() if tabcol.unit is not None else ""
        pfcol = fits.Column(
            name=tabcol.name, unit=unit, array=tabcol.data, format=col_format
        )
        fits_columns.append(pfcol)
        col_descriptions.append(tabcol.description)

    tablehdu = fits.BinTableHDU.from_columns(fits_columns)
    tablehdu.writeto(filename, overwrite=overwrite)


@delayed()
def convert_pdf(pdf_file: str, plots_dir: str, spec_dir: str) -> None:
    """Convert a PDF to a PNG

    Args:
        pdf_file (str): PDF file to convert
    """
    png_file = pdf_file.replace(".pdf", "")
    png_file = os.path.join(plots_dir, os.path.basename(png_file))
    logger.info(f"Converting {pdf_file} to {png_file}")

    cmd = f"pdftoppm {pdf_file} {png_file} -png"
    logger.info(cmd)
    sp.run(cmd.split(), check=True)
    # Grr, pdftoppm doesn't preserve the file name
    actual_name = f"{png_file}-1.png"
    wanted_name = f"{png_file}.png"
    if os.path.exists(actual_name):
        os.rename(actual_name, wanted_name)
    assert os.path.exists(wanted_name)

    # Make thumbnail plot for CASDA - needs same name as FITS
    thumb = os.path.join(spec_dir, os.path.basename(wanted_name))
    rename_dict = {
        "_spectra-plots.png": "_polspec.png",
        # "_RMSF-dirtyFDF-plots.png": "_FDFdirty.png",
        # "_cleanFDF-plots.png": "_FDFclean.png",
    }
    for old, new in rename_dict.items():
        if old in thumb:
            thumb = thumb.replace(old, new)
            try_symlink(os.path.abspath(wanted_name), thumb)
    # other_rename = {
    #     "_spectra.png": "_noise.png",
    #     "_FDFdirty.png": "_RMSF.png",
    #     "_FDFclean.png": "_FDFmodel.png",
    # }
    # for old, new in other_rename.items():
    #     if old in thumb:
    #         thumb = thumb.replace(old, new)
    # try_symlink(os.path.abspath(wanted_name), thumb)


def find_plots(data_dir: str = ".") -> list:
    """Find plots in a directory

    Args:
        data_dir (str, optional): Data containg cutouts directory. Defaults to ".".

    Returns:
        list: List of plots
    """
    cut_dir = os.path.join(data_dir, "cutouts")
    logger.info(f"Globbing for plots in {cut_dir}")
    plots = glob(os.path.join(os.path.join(cut_dir, "RACS_*"), "*.pdf"))
    logger.info(f"Found {len(plots)} plots")
    return plots


def main(
    polcatf: str,
    data_dir: str = ".",
    prep_type: str = "test",
    do_update_cubes: bool = False,
    do_convert_spectra: bool = False,
    do_convert_plots: bool = False,
    verbose: bool = False,
    batch_size: int = 10,
    outdir=None,
):
    """Main function"""
    # Re-register astropy units
    for unit in (u.deg, u.hour, u.hourangle, u.Jy, u.arcsec, u.arcmin, u.beam):
        get_current_unit_registry().add_enabled_units([unit])
    logger.info("Starting")
    logger.info(f"Reading {polcatf}")

    polcat = Table.read(polcatf)
    df = polcat.to_pandas()
    df = df.sort_values(["stokesI_fit_flag", "snr_polint"], ascending=[True, False])
    polcat = polcat[df.index.values]
    polcat.add_index("cat_id")

    logger.info(f"Preparing data for {prep_type} CASDA upload")

    if prep_type == "full":
        pass

    elif prep_type == "cut":
        cut_idx = ~polcat["stokesI_fit_flag"]
        polcat = polcat[cut_idx]

    elif prep_type == "test":
        polcat = polcat[:100]

    else:
        msg = f"Unknown prep_type: {prep_type}"
        raise ValueError(msg)

    casda_dir = os.path.join(data_dir, f"casda_{prep_type}")
    try_mkdir(casda_dir)

    # Link catalgoue to casda directory
    cat_dir = os.path.join(casda_dir, "catalogues")
    try_mkdir(cat_dir)
    if prep_type != "full":
        out_cat = os.path.join(
            cat_dir, os.path.basename(polcatf).replace(".xml", f".{prep_type}.xml")
        )
        write_votable(polcat, out_cat)
    else:
        try_symlink(
            os.path.abspath(polcatf), os.path.join(cat_dir, os.path.basename(polcatf))
        )

    cube_outputs = []
    if do_update_cubes:
        logger.info("Updating cubelets")
        cube_dir = os.path.join(casda_dir, "cubelets")
        try_mkdir(cube_dir)

        cubes = find_cubes(data_dir=data_dir)
        # Check if we have a cube for each source
        try:
            assert len(cubes) == len(
                set(polcat["source_id"])
            ), "Number of cubes does not match number of sources"
        except AssertionError:
            logger.warning(
                f"Found {len(cubes)} cubes, expected {len(set(polcat['source_id']))}"
            )
            if len(cubes) < len(set(polcat["source_id"])):
                logger.critical("Some cubes are missing on disk!")
                raise
            else:
                logger.warning("Need to exclude some cubes")
                source_ids = []
                for i, cube in enumerate(cubes):
                    basename = os.path.basename(cube)
                    cut_idx = basename.find(".cutout.")
                    source_id = basename[:cut_idx]
                    source_ids.append(source_id)
                in_idx = np.isin(source_ids, polcat["source_id"])
                cubes = list(np.array(cubes)[in_idx])
                logger.warning(
                    f"I had to exclude {np.sum(~in_idx)} sources that were not in the catalogue"
                )
                # Write missing source IDs to disk
                rem_ids = list(set(np.array(source_ids)[~in_idx]))
                outf = os.path.join(casda_dir, "excluded_sources.txt")
                logger.info(f"Writing excluded source IDs to {outf}")
                with open(outf, "w") as f:
                    for rid in rem_ids:
                        f.write(f"{rid}\n")
                assert (
                    len(cubes) == len(set(polcat["source_id"]))
                ), f"Number of cubes does not match number of sources -- {len(cubes)=} and {len(set(polcat['source_id']))=}"

        unique_ids, unique_idx = np.unique(polcat["source_id"], return_index=True)
        lookup = dict(zip(unique_ids, unique_idx))
        with tqdm(total=len(cubes), desc="Sorting cubes", file=TQDM_OUT) as pbar:

            def my_sorter(x, lookup=lookup, pbar=pbar):
                basename = x.split("/")[-1]
                cut_idx = basename.find(".cutout")
                sourcename = basename[:cut_idx]
                idx = lookup[sourcename]
                if "weight" in basename:
                    idx += 3
                for i, s in enumerate((".i.", ".q.", ".u.")):
                    if s in basename:
                        idx += i
                        pbar.update(1)
                        return idx
                return None

            sorted_cubes = sorted(cubes, key=my_sorter)

        for cube in sorted_cubes:
            out = update_cube(cube=cube, cube_dir=cube_dir)
            cube_outputs.append(out)

    spectra_outputs = []
    if do_convert_spectra:
        logger.info("Converting spectra")
        spec_dir = os.path.join(casda_dir, "spectra")
        try_mkdir(spec_dir)
        spectra = find_spectra(data_dir=data_dir)
        if prep_type != "full":
            # Drop spectra that are not in the cut catalogue
            cat_ids = []
            for i, spec in enumerate(spectra):
                basename = os.path.basename(spec)
                cut_idx = basename.find(".dat")
                cat_id = basename[:cut_idx]
                cat_ids.append(cat_id)
            in_idx = np.isin(cat_ids, polcat["cat_id"])
            spectra = list(np.array(spectra)[in_idx])
        assert len(spectra) == len(
            polcat
        ), f"{len(spectra)=} and {len(polcat)=}"  # Sanity check

        unique_ids, unique_idx = np.unique(polcat["cat_id"], return_index=True)
        lookup = dict(zip(unique_ids, unique_idx))
        with tqdm(total=len(spectra), desc="Sorting spectra", file=TQDM_OUT) as pbar:

            def my_sorter(x, lookup=lookup, pbar=pbar):
                basename = x.split("/")[-1]
                sourcename = basename.replace(".dat", "")
                idx = lookup[sourcename]
                pbar.update(1)
                return idx

            sorted_spectra = sorted(spectra, key=my_sorter)

        # Loop over spectra and convert to FITS
        gauss_ids = [
            os.path.basename(spectrum).replace(".dat", "")
            for spectrum in sorted_spectra
        ]
        polcat = polcat.loc[gauss_ids]
        # hash filename using current time and hashlib
        fname_polcat_hash = (
            f"polcat_{hashlib.sha256(str(time.time()).encode()).hexdigest()}.pkl"
        )
        with open(fname_polcat_hash, "wb") as f:
            pickle.dump(polcat, f)
        for i, (spectrum, gauss_id, row) in enumerate(
            tqdm(
                zip(sorted_spectra, gauss_ids, polcat),
                total=len(sorted_spectra),
                desc="Converting spectra",
                file=TQDM_OUT,
            )
        ):
            assert gauss_id == row["cat_id"]
            spectrum_table = convert_spectra(
                number=i,
                spectrum=spectrum,
                fname_polcat_hash=fname_polcat_hash,
                spec_dir=spec_dir,
            )
            spectra_outputs.append(spectrum_table)

    plot_outputs = []
    if do_convert_plots:
        logger.info("Converting plots")
        plots_dir = os.path.join(casda_dir, "plots")
        try_mkdir(plots_dir)
        spec_dir = os.path.join(casda_dir, "spectra")
        try_mkdir(spec_dir)
        plots = find_plots(data_dir=data_dir)
        unique_ids, unique_idx = np.unique(polcat["cat_id"], return_index=True)
        lookup = dict(zip(unique_ids, unique_idx))

        if prep_type != "full":
            # Drop plots that are not in the cut catalogue
            cat_ids = []
            for i, plot in enumerate(plots):
                basename = os.path.basename(plot)
                for plot_type in ("_spectra", "_RMSF", "_clean"):
                    if plot_type in basename:
                        cat_id = basename[: basename.find(plot_type)]
                cat_ids.append(cat_id)
            in_idx = np.isin(cat_ids, polcat["cat_id"])
            plots = list(np.array(plots)[in_idx])

        assert (
            len(plots) == len(polcat) * 3
        ), f"{len(plots)=} and {len(polcat)=}"  # Sanity check

        with tqdm(total=len(plots), desc="Sorting plots", file=TQDM_OUT) as pbar:

            def my_sorter(x, lookup=lookup, pbar=pbar):
                fname = x.split("/")[-1]
                for plot_type in ("_spectra", "_RMSF", "_clean"):
                    if plot_type in fname:
                        sourcename = fname[: fname.find(plot_type)]
                        break
                idx = lookup[sourcename]
                pbar.update(1)
                return idx

            sorted_plots = sorted(plots, key=my_sorter)
        for plot in sorted_plots:
            out = convert_pdf(pdf_file=plot, plots_dir=plots_dir, spec_dir=spec_dir)
            plot_outputs.append(out)

    for name, outputs in zip(
        ("cubes", "spectra", "plots"), (cube_outputs, spectra_outputs, plot_outputs)
    ):
        logger.info(f"Starting work on {len(outputs)} {name}")

        futures = chunk_dask(
            outputs=outputs,
            task_name=name,
            progress_text=f"Preparing {name} for CASDA",
            verbose=verbose,
            batch_size=batch_size,
        )
        if name == "spectra" and len(outputs) > 0:
            client = get_dask_client()
            # Get concrete results
            spectrum_tables = client.gather(client.compute(futures))
            # Add all spectrum_tables to a tar ball
            tarball = os.path.join(casda_dir, f"spice_racs_dr1_polspec_{prep_type}.tar")
            logger.info(f"Adding spectra to tarball {tarball}")
            with tarfile.open(tarball, "w") as tar:
                for spectrum_table in tqdm(
                    spectrum_tables,
                    "Adding spectra to tarball",
                    file=TQDM_OUT,
                ):
                    tar.add(spectrum_table, arcname=os.path.basename(spectrum_table))

    if do_convert_spectra:
        os.remove(fname_polcat_hash)

    logger.info("Done")


def cli():
    """Command line interface"""
    parser = argparse.ArgumentParser(description="Prepare files for CASDA upload")
    parser.add_argument(
        "data_dir",
        help="Directory containing the 'cutouts' to be corrected/uploaded",
        type=str,
        default=".",
        metavar="data_dir",
    )
    parser.add_argument(
        "polcat",
        help="Path to the polarisation catalogue.",
        type=str,
        metavar="polcat",
    )

    parser.add_argument(
        "prep_type",
        choices=["full", "cut", "test"],
        help="Type of data to prepare",
        type=str,
        metavar="prep_type",
        default="test",
    )

    parser.add_argument(
        "--convert-cubes", action="store_true", help="Update cubes", default=False
    )
    parser.add_argument(
        "--convert-spectra", action="store_true", help="Convert spectra", default=False
    )
    parser.add_argument(
        "--convert-plots", action="store_true", help="Convert plots", default=False
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Verbose output",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Debug output",
    )
    parser.add_argument(
        "--mpi",
        action="store_true",
        help="Use MPI",
    )
    parser.add_argument(
        "--batch_size",
        help="Number parallel jobs to run",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--interface",
        help="Interface to use",
        type=str,
        default="ipogif0",
    )
    parser.add_argument(
        "--outdir",
        help="Output directory",
        type=str,
        default=None,
    )
    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    main(
        polcatf=args.polcat,
        data_dir=args.data_dir,
        prep_type=args.prep_type,
        do_update_cubes=args.convert_cubes,
        do_convert_spectra=args.convert_spectra,
        do_convert_plots=args.convert_plots,
        verbose=args.verbose,
        batch_size=args.batch_size,
        outdir=args.outdir,
    )


if __name__ == "__main__":
    cli()
