#!/usr/bin/env python3
"""Prepare files for CASDA upload"""
import argparse
import logging as log
import os
import subprocess as sp
import time
import traceback
from glob import glob
from typing import Dict, List, Tuple

import astropy.units as u
import dask.array as da
import dask.bag as db
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polspectra
from astropy.io import fits
from astropy.table import Column, Row, Table, vstack
from astropy.visualization import (
    ImageNormalize,
    LogStretch,
    MinMaxInterval,
    SqrtStretch,
    ZScaleInterval,
)
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from dask import delayed
from dask.delayed import Delayed
from dask.distributed import Client, LocalCluster
from dask_mpi import initialize
from IPython import embed
from radio_beam import Beam
from tqdm.auto import tqdm, trange

from spiceracs.utils import chunk_dask, tqdm_dask, try_mkdir, try_symlink, zip_equal


def make_thumbnail(cube_f: str, cube_dir: str):
    cube = np.squeeze(fits.getdata(cube_f))
    head = fits.getheader(cube_f)
    wcs = WCS(head)
    med = np.nanmedian(cube, axis=0)
    med_wcs = wcs.celestial
    pix_scales = proj_plane_pixel_scales(med_wcs) * u.deg
    beam = Beam.from_fits_header(head)
    ellipse = beam.ellipse_to_plot(
        xcen=10,
        ycen=10,
        pixscale=pix_scales[0],  # Assume same pixel scale in x and y
    )
    fig = plt.figure(facecolor="w")
    if ".weights." in cube_f:
        cmap = "gray"
        norm = ImageNormalize(med, interval=MinMaxInterval())
    else:
        if ".i." in cube_f:
            cmap = "viridis"
            norm = ImageNormalize(med, interval=MinMaxInterval())
        else:
            cmap = "coolwarm"
            absmax = np.max(np.abs(med))
            norm = ImageNormalize(med, vmin=-absmax, vmax=absmax)
    ax = fig.add_subplot(111, projection=med_wcs)
    im = ax.imshow(med, origin="lower", cmap=cmap, norm=norm)
    ax.add_patch(ellipse)
    ellipse.set_clip_box(ax.bbox)
    ellipse.set_facecolor("none")
    ellipse.set_edgecolor("k")
    ellipse.set_hatch("//")
    ax.set_xlabel("RA")
    ax.set_ylabel("Dec")
    fig.colorbar(im, ax=ax, label=f"{u.Unit(head['BUNIT']):latex_inline}")
    outf = os.path.join(cube_dir, os.path.basename(cube_f).replace(".fits", ".png"))
    log.info(f"Saving thumbnail to {outf}")
    fig.savefig(outf, dpi=300)
    plt.close(fig)


def find_spectra(data_dir: str = ".") -> list:
    """Find spectra in from cutouts directory

    Args:
        data_dir (str, optional): Directory containing cutouts directory. Defaults to ".".

    Returns:
        list: List of spectra in ascii format
    """
    cut_dir = os.path.join(data_dir, "cutouts")
    log.info(f"Globbing for spectra in {cut_dir}")
    spectra = glob(os.path.join(os.path.join(cut_dir, "*"), "*[0-9].dat"))
    log.info(f"Found {len(spectra)} spectra (in frequency space)")
    return spectra


@delayed()
def convert_spectra(
    number: int,
    spectrum: str,
    cat_row: Row,
    spec_dir: str = ".",
) -> str:
    """Convert a ascii spectrum to FITS

    Args:
        spectrum (str): Name of ASCII spectrum file
        spec_dir (str, optional): Directory to save FITS spectrum. Defaults to '.'.
    """
    cat_row = Table(cat_row)
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

    spectrum_table = polspectra.from_arrays(
        long_array=cat_row["ra"],
        lat_array=cat_row["dec"],
        freq_array=freq,
        stokesI=(full_data.I.values)[np.newaxis, :],
        stokesI_error=(full_data.dI.values)[np.newaxis, :],
        stokesQ=(full_data.Q.values)[np.newaxis, :],
        stokesQ_error=(full_data.dQ.values)[np.newaxis, :],
        stokesU=(full_data.Q.values)[np.newaxis, :],
        stokesU_error=(full_data.dQ.values)[np.newaxis, :],
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

    for col, desc in tqdm(pol_df_cols.items(), desc="Adding spectrum columns"):
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
            raise ValueError(f"Unknown column {col}")

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
        spec_dir, os.path.basename(spectrum).replace(".dat", f"_polspec.fits")
    )
    log.info(f"Saving to {outf}")
    spectrum_table.write_FITS(outf, overwrite=True)
    # Add object to header
    with fits.open(outf, mode="update") as hdul:
        hdul[0].header["OBJECT"] = (cat_row["cat_id"][0], "Gaussian ID")
        hdul.flush()

    return outf


@delayed()
def update_cube(cube: str, cube_dir: str) -> None:
    """Update cube headers and symlink to CASDA area

    Args:
        cube (str): Cubelet path
        cube_dir (str): CASDA cublet directory
    """
    # Get source ID and update header
    source_id = os.path.basename(os.path.dirname(cube))
    with fits.open(cube, mode="update") as hdul:
        hdul[0].header["OBJECT"] = (source_id, "Source ID")
        hdul.flush()

    # Move cube to cubelets directory
    src = os.path.abspath(cube)
    dst = os.path.join(cube_dir, os.path.basename(cube))
    # Fix old RACS name
    dst = dst.replace("RACS_test4_1.05_", "RACS_")
    log.info(f"Copying {src} to {dst}")
    try_symlink(src, dst)
    make_thumbnail(dst, cube_dir)


def find_cubes(data_dir: str = ".") -> list:
    """Find cubelets in a directory

    Args:
        data_dir (str, optional): Data containg cutouts directory. Defaults to ".".

    Returns:
        list: List of cubelets
    """
    cut_dir = os.path.join(data_dir, "cutouts")
    log.info(f"Globbing for cubes in {cut_dir}")
    cubes = glob(os.path.join(os.path.join(cut_dir, "*"), "*.linmos.fits"))
    log.info(f"Found {len(cubes)} cubes")
    return cubes


@delayed(nout=2)
def init_polspec(
    casda_dir: str,
    spectrum_table_0: str,
    outdir: str = None,
):
    if outdir is None:
        outdir = casda_dir
    polspec_0 = polspectra.from_FITS(spectrum_table_0)
    out_fits = os.path.join(os.path.abspath(outdir), "spice_racs_dr1_polspec.fits")
    log.info(f"Saving to {out_fits}")
    polspec_0.write_FITS(out_fits, overwrite=True)

    out_hdf = os.path.join(os.path.abspath(outdir), "spice_racs_dr1_polspec.hdf5")
    log.info(f"Saving to {out_hdf}")
    polspec_0.write_HDF5(out_hdf, overwrite=True, compress=True)

    return out_fits, out_hdf


# Reworked from polspectra.py
def write_polspec(table: Table, filename: str, overwrite: bool = False):
    """Write the polspectra to a FITS file.
    Parameters:
        table: astropy.table.Table
        filename : str
            Name and relative path of the file to save to.
        overwrite : bool [False]
            Overwrite the file if it already exists?"""
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
        if tabcol.unit != None:
            unit = tabcol.unit.to_string()
        else:
            unit = ""
        pfcol = fits.Column(
            name=tabcol.name, unit=unit, array=tabcol.data, format=col_format
        )
        fits_columns.append(pfcol)
        col_descriptions.append(tabcol.description)

    tablehdu = fits.BinTableHDU.from_columns(fits_columns)
    tablehdu.writeto(filename, overwrite=overwrite)


@delayed()
def add_polspec_row(
    out_fits: str,
    out_hdf: str,
    spectrum_table: polspectra.polarizationspectra,
) -> None:
    # Add row to FITS table
    with fits.open(out_fits, mode="denywrite", memmap=True) as hdul:
        data = hdul[1].data
        header = hdul[1].header
        table = Table(data)
        for i in range(1, header["TFIELDS"] + 1):
            if "TUNIT{}".format(i) in header.keys():
                table[header["TTYPE{}".format(i)]].unit = header["TUNIT{}".format(i)]
    stack = vstack([table, spectrum_table.table], join_type="exact")
    write_polspec(stack, out_fits, overwrite=True)

    # TODO: Make this work
    # # Add row to HDF5 table
    # stack = polspectra.polarizationspectra()
    # stack.read_HDF5(out_hdf)
    # stack.merge_tables(spectrum_table)
    # stack.write_HDF5(out_hdf, overwrite=True, compress=True)


@delayed()
def convert_pdf(pdf_file: str, plots_dir: str, spec_dir: str) -> None:
    """Convert a PDF to a PNG

    Args:
        pdf_file (str): PDF file to convert
    """
    png_file = pdf_file.replace(".pdf", "")
    png_file = os.path.join(plots_dir, os.path.basename(png_file))
    log.info(f"Converting {pdf_file} to {png_file}")

    cmd = f"pdftoppm {pdf_file} {png_file} -png"
    log.info(cmd)
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
        "_spectra-plots.png": "_spectra.png",
        "_RMSF-dirtyFDF-plots.png": "_FDFdirty.png",
        "_cleanFDF-plots.png": "_FDFclean.png",
    }
    for old, new in rename_dict.items():
        if old in thumb:
            thumb = thumb.replace(old, new)
    try_symlink(os.path.abspath(wanted_name), thumb)
    other_rename = {
        "_spectra.png": "_noise.png",
        "_FDFdirty.png": "_RMSF.png",
        "_FDFclean.png": "_FDFmodel.png",
    }
    for old, new in other_rename.items():
        if old in thumb:
            thumb = thumb.replace(old, new)
    try_symlink(os.path.abspath(wanted_name), thumb)


def find_plots(data_dir: str = ".") -> list:
    """Find plots in a directory

    Args:
        data_dir (str, optional): Data containg cutouts directory. Defaults to ".".

    Returns:
        list: List of plots
    """
    cut_dir = os.path.join(data_dir, "cutouts")
    log.info(f"Globbing for plots in {cut_dir}")
    plots = glob(os.path.join(os.path.join(cut_dir, "RACS_*"), "*.pdf"))
    log.info(f"Found {len(plots)} plots")
    return plots


def main(
    polcatf: str,
    client: Client,
    data_dir: str = ".",
    do_update_cubes: bool = False,
    do_convert_spectra: bool = False,
    do_convert_plots: bool = False,
    verbose: bool = False,
    test: bool = False,
    batch_size: int = 10,
    outdir=None,
):
    """Main function"""

    log.info("Starting")
    log.info(f"Dask client: {client}")
    log.info(f"Reading {polcatf}")
    polcat = Table.read(polcatf)
    df = polcat.to_pandas()
    df = df.sort_values(["stokes_I_fit_flag", "snr_polint"], ascending=[True, False])
    polcat = polcat[df.index.values]
    polcat.add_index("cat_id")

    casda_dir = (
        os.path.join(data_dir, "casda")
        if not test
        else os.path.join(data_dir, "casda_test")
    )
    try_mkdir(casda_dir)

    # Link catalgoue to casda directory
    cat_dir = os.path.join(casda_dir, "catalogues")
    try_mkdir(cat_dir)
    try_symlink(
        os.path.abspath(polcatf), os.path.join(cat_dir, os.path.basename(polcatf))
    )

    cube_outputs = []
    if do_update_cubes:
        log.info("Updating cubelets")
        cube_dir = os.path.join(casda_dir, "cubelets")
        try_mkdir(cube_dir)

        cubes = find_cubes(data_dir=data_dir)
        # Check if we have a cube for each source
        try:
            assert (
                len(cubes) == len(set(polcat["source_id"])) * 3 * 2
            )  # 3 pols x (image, weight)
        except AssertionError:
            log.warning(
                f"Found {len(cubes)} cubes, expected {len(set(polcat['source_id'])) * 3 * 2}"
            )
            if len(cubes) < len(set(polcat["source_id"])) * 3 * 2:
                log.critical("Some cubes are missing on disk!")
                raise
            else:
                log.warning("Need to exclude some cubes")
                source_ids = []
                for i, cube in enumerate(cubes):
                    basename = os.path.basename(cube)
                    cut_idx = basename.find(".cutout.")
                    source_id = basename[:cut_idx]
                    source_ids.append(source_id)
                in_idx = np.isin(source_ids, polcat["source_id"])
                cubes = list(np.array(cubes)[in_idx])
                log.warning(
                    f"I had to exclude {np.sum(~in_idx)} sources that were not in the catalogue"
                )
                # Write missing source IDs to disk
                rem_ids = list(set(np.array(source_ids)[~in_idx]))
                outf = os.path.join(casda_dir, "excluded_sources.txt")
                log.info(f"Writing excluded source IDs to {outf}")
                with open(outf, "w") as f:
                    for rid in rem_ids:
                        f.write(f"{rid}\n")
                assert len(cubes) == len(set(polcat["source_id"])) * 3 * 2

        unique_ids, unique_idx = np.unique(polcat["source_id"], return_index=True)
        lookup = {sid: i for sid, i in zip(unique_ids, unique_idx)}
        with tqdm(total=len(cubes), desc="Sorting cubes") as pbar:

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

            sorted_cubes = sorted(cubes, key=my_sorter)

        for cube in sorted_cubes:
            out = update_cube(cube=cube, cube_dir=cube_dir)
            cube_outputs.append(out)

    spectrum_tables = []
    spectra_outputs = []
    if do_convert_spectra:
        log.info("Converting spectra")
        spec_dir = os.path.join(casda_dir, "spectra")
        try_mkdir(spec_dir)
        spectra = find_spectra(data_dir=data_dir)
        assert len(spectra) == len(polcat)  # Sanity check

        unique_ids, unique_idx = np.unique(polcat["cat_id"], return_index=True)
        lookup = {sid: i for sid, i in zip(unique_ids, unique_idx)}
        with tqdm(total=len(spectra), desc="Sorting spectra") as pbar:

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
        sorted_polcat = polcat.loc[gauss_ids]
        for i, (spectrum, gauss_id, row) in enumerate(
            tqdm(
                zip(sorted_spectra, gauss_ids, sorted_polcat),
                total=len(sorted_spectra),
                desc="Converting spectra",
            )
        ):
            assert gauss_id == row["cat_id"]
            spectrum_table = convert_spectra(
                number=i,
                spectrum=spectrum,
                cat_row=row,
                spec_dir=spec_dir,
            )
            spectrum_tables.append(spectrum_table)
        # Init spectrum table
        out_fits, out_hdf = init_polspec(
            casda_dir=casda_dir,
            spectrum_table_0=spectrum_tables[0],
            outdir=outdir,
        )

        for spectrum_table in spectrum_tables:
            out = add_polspec_row(
                out_fits=out_fits,
                out_hdf=out_hdf,
                spectrum_table=spectrum_table,
            )
            spectra_outputs.append(out)

    plot_outputs = []
    if do_convert_plots:
        log.info("Converting plots")
        plots_dir = os.path.join(casda_dir, "plots")
        try_mkdir(plots_dir)
        spec_dir = os.path.join(casda_dir, "spectra")
        try_mkdir(spec_dir)
        plots = find_plots(data_dir=data_dir)
        unique_ids, unique_idx = np.unique(polcat["cat_id"], return_index=True)
        lookup = {sid: i for sid, i in zip(unique_ids, unique_idx)}
        with tqdm(total=len(plots), desc="Sorting plots") as pbar:

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
        if test:
            if name == "spectra":
                n_things = 10 * 10
            elif name == "cubes":
                n_things = 60 * 10
            else:
                n_things = 30 * 10
            outputs = outputs[:n_things]
            log.info(f"Testing {len(outputs)} {name}")
        else:
            log.info(f"Starting work on {len(outputs)} {name}")

        futures = chunk_dask(
            outputs=outputs,
            client=client,
            task_name=name,
            progress_text=f"Preparing {name} for CASDA",
            verbose=verbose,
            batch_size=batch_size,
        )

    log.info("Done")


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
        "--update-cubes", action="store_true", help="Update cubes", default=False
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
        "--test",
        action="store_true",
        help="Test mode",
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
        log.basicConfig(
            level=log.INFO,
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    elif args.debug:
        log.basicConfig(
            level=log.DEBUG,
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    else:
        log.basicConfig(
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

    if args.mpi:
        initialize(
            interface=args.interface,
            local_directory="/dev/shm",
        )
        cluster = None
    else:
        cluster = LocalCluster(
            n_workers=12,
            processes=True,
            threads_per_worker=1,
            local_directory="/dev/shm",
        )
    with Client(cluster) as client:
        log.debug(client)
        main(
            polcatf=args.polcat,
            client=client,
            data_dir=args.data_dir,
            do_update_cubes=args.update_cubes,
            do_convert_spectra=args.convert_spectra,
            do_convert_plots=args.convert_plots,
            verbose=args.verbose,
            test=args.test,
            batch_size=args.batch_size,
            outdir=args.outdir,
        )


if __name__ == "__main__":
    cli()
