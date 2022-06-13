#!/usr/bin/env python3
"""Prepare files for CASDA upload"""
import os
import time
from typing import Tuple, List, Dict
from spiceracs.utils import try_mkdir, tqdm_dask, try_symlink, chunk_dask
import polspectra
from glob import glob
from astropy.io import fits
from astropy.table import Table, Column
import astropy.units as u
import numpy as np
import logging as log
import traceback
import argparse
from dask import delayed
from dask.distributed import Client, LocalCluster
import dask.bag as db
import dask.array as da
import subprocess as sp
from dask_mpi import initialize
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize, ZScaleInterval, LogStretch)
import astropy.units as u
from radio_beam import Beam
from tqdm.auto import tqdm, trange
import pandas as pd
from IPython import embed

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
        pixscale=pix_scales[0], # Assume same pixel scale in x and y
        
    )
    fig = plt.figure(facecolor='w')
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
    im = ax.imshow(med, origin='lower', cmap=cmap, norm=norm)
    ax.add_patch(ellipse)
    ellipse.set_clip_box(ax.bbox)
    ellipse.set_facecolor('none')
    ellipse.set_edgecolor('k')
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


@delayed
def convert_spectra(
    spectrum: str, 
    ra: float,
    dec: float,
    gauss_id: str,
    spec_dir: str = ".",
    ) -> dict:
    """Convert a ascii spectrum to FITS

    Args:
        spectrum (str): Name of ASCII spectrum file
        spec_dir (str, optional): Directory to save FITS spectrum. Defaults to '.'.
    """
    rmsf = u.def_unit("RMSF")

    # First deal with the frequency data
    full_data = np.loadtxt(spectrum).T
    freq, i_data, q_data, u_data, i_noise, q_noise, u_noise = full_data

    freq = freq * u.Hz

    # Hard code the pixel size for now
    pix_size = 2.5 * u.arcsec
    log.warning(f"Using a pixel size of {pix_size}")

    data = np.array([i_data, q_data, u_data]) * u.Jy / u.beam
    noise = np.array([i_noise, q_noise, u_noise]) * u.Jy / u.beam
    # Add dummy axes for RA/DEC
    data = data[:, :, np.newaxis, np.newaxis]
    noise = noise[:, :, np.newaxis, np.newaxis]

    # Free up memory
    del i_data, q_data, u_data, i_noise, q_noise, u_noise, full_data
    # Create the data FITS files
    for d, n in zip((data, noise), ("spectra", "noise")):
        hdu = fits.PrimaryHDU(data=data.value)
        hdu.header["NAXIS"] = len(d.shape)
        hdu.header["NAXIS1"] = 1
        hdu.header["NAXIS2"] = 1
        hdu.header["NAXIS3"] = len(freq)
        hdu.header["NAXIS4"] = len(d)
        hdu.header["CTYPE1"] = "RA---SIN"
        hdu.header["CTYPE2"] = "DEC--SIN"
        hdu.header["CTYPE3"] = "FREQ"
        hdu.header["CTYPE4"] = "STOKES"
        hdu.header["CRVAL1"] = ra
        hdu.header["CRVAL2"] = dec
        hdu.header["CRVAL3"] = freq[0].value
        hdu.header["CRVAL4"] = 1
        hdu.header["CRPIX1"] = 1
        hdu.header["CRPIX2"] = 1
        hdu.header["CRPIX3"] = 1
        hdu.header["CRPIX4"] = 1
        hdu.header["CDELT1"] = pix_size.to(u.deg).value
        hdu.header["CDELT2"] = pix_size.to(u.deg).value
        hdu.header["CDELT3"] = np.diff(freq)[0].value
        hdu.header["CDELT4"] = 1
        hdu.header["CUNIT1"] = u.deg.to_string(format="fits")
        hdu.header["CUNIT2"] = u.deg.to_string(format="fits")
        hdu.header["CUNIT3"] = freq.unit.to_string(format="fits")
        hdu.header["CUNIT4"] = "STOKES"
        hdu.header["BUNIT"] = d.unit.to_string(format="fits")
        hdu.header["OBJECT"] = (gauss_id, "Gaussian ID")

        f = os.path.join(
            spec_dir, os.path.basename(spectrum).replace(".dat", f"_{n}.fits")
        )
        log.info(f"Writing {f} from {spectrum}")
        hdu.writeto(f, overwrite=True)

    # Now deal with the RM data
    suffixes = (
        "_FDFclean.dat",
        "_FDFdirty.dat",
        "_FDFmodel.dat",
        "_RMSF.dat",
    )
    fdf_clean_data = []
    fdf_dirty_data = []
    fdf_model_data = []
    rmsf_data = []
    phis_list = []
    phis_long_list = []

    names = {} # type: Dict[str, str]
    for suffix in suffixes:
        unit = u.Jy / u.beam / rmsf

        is_rmsf = "RMSF" in suffix
        is_clean = "clean" in suffix
        is_dirty = "dirty" in suffix
        is_model = "model" in suffix
        if suffix == "_RMSF.dat":
            unit = u.dimensionless_unscaled


        rm_file = spectrum.replace(".dat", suffix)
        full_rm_data = np.loadtxt(rm_file).T
        phis = full_rm_data[0] * u.rad / u.m ** 2
        if is_rmsf:
            phis_long_list.append(phis)
        elif is_clean:
            phis_list.append(phis)
        rm_q_data = full_rm_data[1]
        rm_u_data = full_rm_data[2]
        rm_data = np.array([rm_q_data, rm_u_data]) * unit
        # Add dummy axis for RA/DEC
        rm_data = rm_data[:, :, np.newaxis, np.newaxis]
        if is_rmsf:
            rmsf_data.append(rm_data)
        elif is_clean:
            fdf_clean_data.append(rm_data)
        elif is_dirty:
            fdf_dirty_data.append(rm_data)
        elif is_model:
            fdf_model_data.append(rm_data)

        # Create HDU
        hdu = fits.PrimaryHDU(data=rm_data.value)
        hdu.header["NAXIS"] = len(rm_data.shape)
        hdu.header["NAXIS1"] = 1
        hdu.header["NAXIS2"] = 1
        hdu.header["NAXIS3"] = len(phis)
        hdu.header["NAXIS4"] = len(rm_data)
        hdu.header["CTYPE1"] = "RA---SIN"
        hdu.header["CTYPE2"] = "DEC--SIN"
        hdu.header["CTYPE3"] = "FDEP"
        hdu.header["CTYPE4"] = "STOKES"
        hdu.header["CRVAL1"] = ra
        hdu.header["CRVAL2"] = dec
        hdu.header["CRVAL3"] = phis[0].value
        hdu.header["CRVAL4"] = 2  # Stokes Q
        hdu.header["CRPIX1"] = 1
        hdu.header["CRPIX2"] = 1
        hdu.header["CRPIX3"] = 1
        hdu.header["CRPIX4"] = 1
        hdu.header["CDELT1"] = pix_size.to(u.deg).value
        hdu.header["CDELT2"] = pix_size.to(u.deg).value
        hdu.header["CDELT3"] = np.diff(phis)[0].value
        hdu.header["CDELT4"] = 1
        hdu.header["CUNIT1"] = u.deg.to_string(format="fits")
        hdu.header["CUNIT2"] = u.deg.to_string(format="fits")
        hdu.header["CUNIT3"] = phis.unit.to_string()
        hdu.header["CUNIT4"] = "STOKES"
        hdu.header["BUNIT"] = rm_data.unit.to_string()
        hdu.header["OBJECT"] = (gauss_id, "Gaussian ID")
        rm_f = os.path.join(
            spec_dir, os.path.basename(rm_file).replace(".dat", ".fits")
        )
        log.info(f"Writing {rm_f} from {rm_file}")
        hdu.writeto(rm_f, overwrite=True)
        names[suffix.replace("_","").replace(".dat","")] = rm_f

    freq = np.squeeze(freq)
    data = np.squeeze(data)
    noise = np.squeeze(noise)
    phis = np.squeeze(phis_list)
    phis_long = np.squeeze(phis_long_list)
    fdf_clean_data = np.squeeze(fdf_clean_data)
    fdf_dirty_data = np.squeeze(fdf_dirty_data)
    fdf_model_data = np.squeeze(fdf_model_data)
    rmsf_data = np.squeeze(rmsf_data)
    gauss_id = np.squeeze(gauss_id)

    return dict(
            freq_array=freq, 
            StokesI=data[0],
            StokesQ=data[1],
            StokesU=data[2],
            StokesI_error=noise[0],
            StokesQ_error=noise[1],
            StokesU_error=noise[2],
            faraday_depth=phis, 
            faraday_depth_long=phis_long,
            FDF_Q_clean=fdf_clean_data[0],
            FDF_U_clean=fdf_clean_data[1],
            FDF_Q_dirty=fdf_dirty_data[0],
            FDF_U_dirty=fdf_dirty_data[1],
            FDF_Q_model=fdf_model_data[0],
            FDF_U_model=fdf_model_data[1],
            RMSF_Q=rmsf_data[0],
            RMSF_U=rmsf_data[1],
            cat_id=np.array(gauss_id)
        )


@delayed
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

def make_polspec(
    casda_dir: str,
    polcat: Table,
    pol_df: pd.DataFrame,
    outdir: str=None,
) -> None:
    """Make a PolSpectra table

    Args:
        casda_dir (str): CASDA directory
        polcat (Table): Polarisation catalogue
        freqs (np.ndarray): Array of frequency arrays
        data (np.ndarray): Array of data arrays
        noises (np.ndarray): Array of noise arrays
        gauss_ids (np.ndarray): Array of Gaussian IDs
    """
    polcat.add_index("cat_id")
    freq = pol_df["freq_array"][0]

    # Sort everying by gauss_ids
    polcat = polcat.loc[pol_df["cat_id"].values]
    assert np.array_equal(polcat["cat_id"], pol_df["cat_id"].values)

    spectrum_table = polspectra.from_arrays(
        long_array=polcat["ra"],
        lat_array=polcat["dec"],
        freq_array=freq,
        StokesI=pol_df["StokesI"],
        StokesI_error=pol_df["StokesI_error"],
        StokesQ=pol_df["StokesQ"],
        StokesQ_error=pol_df["StokesQ_error"],
        StokesU=pol_df["StokesU"],
        StokesU_error=pol_df["StokesU_error"],
        source_number_array=range(len(pol_df)),
        cat_id=polcat["cat_id"],
        beam_major=polcat["beam_maj"],
        beam_minor=polcat["beam_min"],
        beam_pa=polcat["beam_pa"],
        coordinate_system="icrs",
        channel_width=np.diff(freq)[0],
        telescope=polcat["telescope"],
        epoch=polcat["epoch"],
        integration_time=polcat["int_time"],
        leakage=polcat["leakage"],
        flux_type=polcat["flux_type"],
    )

    rmsf_unit = u.def_unit("RMSF")
    unit = u.Jy / u.beam / rmsf_unit
    radms = u.radian / u.m**2

    for name, unit, description in zip(
        (
            "faraday_depth", 
            "faraday_depth_long", 
            "FDF_Q_dirty", 
            "FDF_U_dirty",
            "FDF_Q_clean",
            "FDF_U_clean",
            "FDF_Q_model",
            "FDF_U_model",
            "RMSF_Q",
            "RMSF_U",
        ),
        (
            radms,
            radms,
            unit,
            unit,
            unit,
            unit,
            unit,
            unit,
            unit,
            unit,
        ),
        (
            "Faraday depth",
            "Faraday depth (long)",
            "Dirty Stokes Q FDF",
            "Dirty Stokes U FDF",
            "Clean Stokes Q FDF",
            "Clean Stokes U FDF",
            "Model Stokes Q FDF",
            "Model Stokes U FDF",
            "Stokes Q RMSF",
            "Stokes U RMSF",
        )
    ):
        data = pol_df[name]
        new_col = Column(
            name=name,
            unit=unit,
            dtype='object',
            shape=(),
            length=spectrum_table.Nrows,
            description=description,
        )
        new_col[:] = [x for x in data]
        spectrum_table.table.add_column(new_col)

    if outdir is None:
        outdir = casda_dir

    outf = os.path.join(os.path.abspath(outdir), "spice_racs_dr1_polspec.fits",)
    log.info(f"Writing {outf}")
    spectrum_table.write_FITS(outf, overwrite=True)

@delayed
def convert_pdf(pdf_file: str, plots_dir:str, spec_dir: str) -> None:
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
    try_symlink(
        os.path.abspath(wanted_name),
        thumb
    )
    other_rename = {
        "_spectra.png": "_noise.png",
        "_FDFdirty.png": "_RMSF.png",
        "_FDFclean.png": "_FDFmodel.png",
    }
    for old, new in other_rename.items():
        if old in thumb:
            thumb = thumb.replace(old, new)
    try_symlink(
        os.path.abspath(wanted_name),
        thumb
    )


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
    df = df.sort_values(["stokes_I_fit_flag", "snr_polint"], ascending=[True,False])
    polcat = polcat[df.index.values]
    polcat.add_index("cat_id")

    casda_dir = os.path.join(data_dir, "casda")
    try_mkdir(casda_dir)

    # Link catalgoue to casda directory
    cat_dir = os.path.join(casda_dir, "catalogues")
    try_mkdir(cat_dir)
    try_symlink(os.path.abspath(polcatf), os.path.join(cat_dir, os.path.basename(polcatf)))

    cube_outputs = []
    if do_update_cubes:
        log.info("Updating cubelets")
        cube_dir = os.path.join(casda_dir, "cubelets")
        try_mkdir(cube_dir)

        cubes = find_cubes(data_dir=data_dir)
        # Check if we have a cube for each source
        try:
            assert len(cubes) == len(set(polcat["source_id"])) * 3 * 2 # 3 pols x (image, weight)
        except AssertionError:
            log.warning(f"Found {len(cubes)} cubes, expected {len(set(polcat['source_id'])) * 3 * 2}")
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
                log.warning(f"I had to exclude {np.sum(~in_idx)} sources that were not in the catalogue")
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

    spectra_outputs = []
    if do_convert_spectra:
        log.info("Converting spectra")
        spec_dir = os.path.join(casda_dir, "spectra")
        try_mkdir(spec_dir)
        spectra = find_spectra(data_dir=data_dir)
        assert len(spectra) == len(polcat) # Sanity check

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
        gauss_ids = [os.path.basename(spectrum).replace(".dat", "") for spectrum in sorted_spectra]
        sorted_polcat = polcat.loc[gauss_ids]
        for spectrum, gauss_id, row in zip(
                sorted_spectra, gauss_ids, sorted_polcat
        ):
            out = convert_spectra(
                spectrum=spectrum, 
                ra=row["ra"],
                dec=row["dec"],
                gauss_id=gauss_id,
                spec_dir=spec_dir
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
                        sourcename = fname[:fname.find(plot_type)]
                        break
                idx = lookup[sourcename]
                pbar.update(1)
                return idx

            sorted_plots = sorted(plots, key=my_sorter)
        for plot in sorted_plots:
            out = convert_pdf(pdf_file=plot, plots_dir=plots_dir, spec_dir=spec_dir)
            plot_outputs.append(out)

    for name, outputs in zip(("cubes", "spectra", "plots"), (cube_outputs, spectra_outputs, plot_outputs)):
        if test:
            if name == "spectra":
                n_things = 10
            elif name == "cubes":
                n_things = 60
            else:
                n_things = 30
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

        # For spectra, we also want to make a polspec catalogue
        if name == "spectra" and len(outputs) > 0:
            log.info("Making polspec catalogue")
            results = [f.compute() for f in futures]
            pol_df = pd.DataFrame(results)
            pol_df.set_index("cat_id", inplace=True, drop=False)

            # keys = futures[0].compute().keys()
            # out_data_lists = {
            #     key: [] for key in keys
            # } # type: Dict[str, list]

            # for future in futures:
            #     result = future.compute()
            #     for key in keys:
            #         out_data_lists[key].append(result[key])
            
            # out_data_arrs = {} # type: Dict[str, np.ndarray]
            # for key in keys:
            #     out_data_arrs[f"{key}_arr"] = np.squeeze(np.array(out_data_lists[key]))
            # Make polspec catalogue
            make_polspec(
                casda_dir=casda_dir,
                polcat=polcat,
                pol_df=pol_df,
                outdir=outdir,
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
        "-v", "--verbose", action="store_true", help="Verbose output",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Debug output",
    )
    parser.add_argument(
        "--test", action="store_true", help="Test mode",
    )
    parser.add_argument(
        "--mpi", action="store_true", help="Use MPI",
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
        client = Client()
    else:
        cluster = LocalCluster(
            n_workers=12, processes=True, threads_per_worker=1, local_directory="/dev/shm"
        )
        client = Client(cluster)
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
    client.close()

if __name__ == "__main__":
    cli()
