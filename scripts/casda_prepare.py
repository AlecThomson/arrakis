#!/usr/bin/env python3
"""Prepare files for CASDA upload"""
import os
import time
from typing import Tuple, List
from spiceracs.utils import try_mkdir, tqdm_dask, try_symlink
import polspectra
from glob import glob
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import numpy as np
import logging as log
import traceback
import argparse
from dask import delayed
from dask.distributed import Client, LocalCluster
import dask.bag as db
import dask.array as da


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
def convert_spectra(spectrum: str, polcat: Table, spec_dir: str = ".") -> Tuple[np.ndarray]:
    """Convert a ascii spectrum to FITS

    Args:
        spectrum (str): Name of ASCII spectrum file
        spec_dir (str, optional): Directory to save FITS spectrum. Defaults to '.'.
    """
    rmsf = u.def_unit("RMSF")
    gauss_id = os.path.basename(spectrum).replace(".dat", "")
    row = polcat.loc[gauss_id]

    # First deal with the frequency data
    full_data = np.loadtxt(spectrum).T
    freq, i_data, q_data, u_data, i_noise, q_noise, u_noise = full_data

    freq = freq * u.Hz

    data = np.array([i_data, q_data, u_data]) * u.Jy / u.beam
    noise = np.array([i_noise, q_noise, u_noise]) * u.Jy / u.beam
    data = data[:, :, np.newaxis, np.newaxis]
    noise = noise[:, :, np.newaxis, np.newaxis]

    # Create the data FITS file
    data_hdu = fits.PrimaryHDU(data=data.value)
    noise_hdu = fits.PrimaryHDU(data=noise.value)
    # TODO: Add location of source
    for hdu in [data_hdu, noise_hdu]:
        hdu.header["NAXIS"] = len(data.shape)
        hdu.header["NAXIS1"] = 1
        hdu.header["NAXIS2"] = 1
        hdu.header["NAXIS3"] = len(freq)
        hdu.header["NAXIS4"] = len(data)
        hdu.header["CTYPE1"] = "RA---SIN"
        hdu.header["CTYPE2"] = "DEC--SIN"
        hdu.header["CTYPE3"] = "FREQ"
        hdu.header["CTYPE4"] = "STOKES"
        hdu.header["CRVAL1"] = row["ra"].to(u.deg).value
        hdu.header["CRVAL2"] = row["dec"].to(u.deg).value
        hdu.header["CRVAL3"] = freq[0].value
        hdu.header["CRVAL4"] = 1
        hdu.header["CRPIX1"] = 1
        hdu.header["CRPIX2"] = 1
        hdu.header["CRPIX3"] = 1
        hdu.header["CRPIX4"] = 1
        hdu.header["CDELT3"] = np.diff(freq)[0].value
        hdu.header["CDELT4"] = 1
        hdu.header["CUNIT1"] = u.deg.to_string(format="fits")
        hdu.header["CUNIT2"] = u.deg.to_string(format="fits")
        hdu.header["CUNIT3"] = freq.unit.to_string(format="fits")
        hdu.header["CUNIT4"] = "STOKES"
        hdu.header["BUNIT"] = data.unit.to_string(format="fits")
        hdu.header["OBJECT"] = (gauss_id, "Gaussian ID")

    data_f = os.path.join(
        spec_dir, os.path.basename(spectrum).replace(".dat", "_spectra.fits")
    )
    noise_f = os.path.join(
        spec_dir, os.path.basename(spectrum).replace(".dat", "_noise.fits")
    )
    for hdu, f in zip((data_hdu, noise_hdu), (data_f, noise_f)):
        log.info(f"Writing {f} from {spectrum}")
        hdu.writeto(f, overwrite=True)

    # Now deal with the RM data
    suffixes = (
        "_FDFclean.dat",
        "_FDFdirty.dat",
        "_FDFmodel.dat",
        "_RMSF.dat",
    )
    for suffix in suffixes:
        unit = u.Jy / u.beam / rmsf

        if suffix == "_RMSF.dat":
            unit = u.dimensionless_unscaled

        rm_file = spectrum.replace(".dat", suffix)
        full_rm_data = np.loadtxt(rm_file).T
        phis = full_rm_data[0] * u.rad / u.m ** 2
        rm_q_data = full_rm_data[1]
        rm_u_data = full_rm_data[2]
        rm_data = np.array([rm_q_data, rm_u_data]) * unit
        rm_data = rm_data[:, :, np.newaxis, np.newaxis]
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
        hdu.header["CRVAL1"] = row["ra"].to(u.deg).value
        hdu.header["CRVAL2"] = row["dec"].to(u.deg).value
        hdu.header["CRVAL3"] = phis[0].value
        hdu.header["CRVAL4"] = 2  # Stokes Q
        hdu.header["CRPIX1"] = 1
        hdu.header["CRPIX2"] = 1
        hdu.header["CRPIX3"] = 1
        hdu.header["CRPIX4"] = 1
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

    return (freq, data, noise, np.array(gauss_id))


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


def find_cubes(data_dir: str = ".") -> list:
    """Find cubelets in a directory

    Args:
        data_dir (str, optional): Data containg cutouts directory. Defaults to ".".

    Returns:
        list: List of cubelets
    """
    cut_dir = os.path.join(data_dir, "cutouts")
    log.info(f"Globbing for cubes in {cut_dir}")
    cubes = glob(os.path.join(os.path.join(cut_dir, "*"), "*linmos.edge.linmos.fits"))
    log.info(f"Found {len(cubes)} cubes")
    return cubes


@delayed
def make_polspec(
    casda_dir: str,
    polcat: Table,
    freqs: list,
    data: list,
    noises: list,
    gauss_ids: list,
) -> None:
    """Make a PolSpectra table

    Args:
        casda_dir (str): CASDA directory
        polcat (Table): Polarisation catalogue
        freqs (list): list of frequency arrays
        data (list): list of data arrays
        noises (list): list of noise arrays
        gauss_ids (list): list of Gaussian IDs
    """
    freq = freqs[0]

    # Sort everying by gauss_ids
    sort_idx = np.argsort(gauss_ids)
    polcat = polcat[sort_idx]
    data = data[sort_idx]
    noises = noises[sort_idx]
    gauss_ids = gauss_ids[sort_idx]
    assert np.array_equal(polcat["cat_id"], gauss_ids)

    spectrum_table = polspectra.from_arrays(
        long_array=polcat["ra"],
        lat_array=polcat["dec"],
        freq_array=freq,
        StokesI=np.array(data)[:, 0],
        StokesI_error=np.array(noises)[:, 0],
        StokesQ=np.array(data)[:, 1],
        StokesQ_error=np.array(noises)[:, 1],
        StokesU=np.array(data)[:, 2],
        StokesU_error=np.array(noises)[:, 2],
        source_number_array=range(len(data)),
        cat_id=polcat["cat_id"],
        unit=(u.Jy / u.beam).to_string(),
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
    outf = os.path.join(casda_dir, "spice_racs_dr1_polspec.fits",)
    log.info(f"Writing {outf}")
    spectrum_table.write_FITS(outf, overwrite=True)


def main(
    polcat: str,
    client: Client,
    data_dir: str = ".",
    do_update_cubes: bool = False,
    do_convert_spectra: bool = False,
    verbose: bool = False,
):
    """Main function"""

    log.info("Starting")
    log.info(f"Dask client: {client}")

    log.info(f"Reading {polcat}")
    polcat = Table.read(polcat)
    polcat.add_index("cat_id")

    casda_dir = os.path.join(data_dir, "casda")
    try_mkdir(casda_dir)

    cube_outputs = []
    if do_update_cubes:
        log.info("Updating cubelets")
        cube_dir = os.path.join(casda_dir, "cubelets")
        try_mkdir(cube_dir)

        cubes = find_cubes(data_dir=data_dir)
        for cube in cubes:
            out = update_cube(cube=cube, cube_dir=cube_dir)
            cube_outputs.append(out)

    spectra_outputs = []
    if do_convert_spectra:
        log.info("Converting spectra")
        spec_dir = os.path.join(casda_dir, "spectra")
        try_mkdir(spec_dir)
        spectra = find_spectra(data_dir=data_dir)
        # Loop over spectra and convert to FITS
        for spectrum in spectra:
            # freq, datum, noise, gauss_id = convert_spectra(spectrum, spec_dir=spec_dir)
            out = convert_spectra(
                spectrum=spectrum, 
                polcat=polcat, 
                spec_dir=spec_dir
            )
            spectra_outputs.append(out)
            # freqs.append(freq)
            # data.append(datum)
            # noises.append(noise)
            # gauss_ids.append(gauss_id)

        # Convert to dask arrays
        # log.info("Converting to dask arrays")
        # freqs_arr = delayed_to_da(freqs, chunk=1000)
        # data_arr = delayed_to_da(data, chunk=1000)
        # noises_arr = delayed_to_da(noises, chunk=1000)
        # gauss_ids_arr = delayed_to_da(gauss_ids, chunk=1000)

        # # freqs, data, noises, gauss_ids = dummy(
        # #     freqs,
        # #     data,
        # #     noises,
        # #     gauss_ids
        # # )
        # # Make polspec
        # out = make_polspec(
        #     casda_dir=casda_dir,
        #     polcat=polcat,
        #     freqs=freqs_arr,
        #     data=data_arr,
        #     noises=noises_arr,
        #     gauss_ids=gauss_ids_arr,
        # )
        # outputs.append(out)

    from IPython import embed

    embed()
    exit()

    log.info("Starting computation")
    futures = client.persist(outputs)
    # dumb solution for https://github.com/dask/distributed/issues/4831
    log.info("I sleep")
    time.sleep(10)
    log.info("Awake!")
    tqdm_dask(futures, desc="Preparing for CASDA", disable=(not verbose))
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
        "-v", "--verbose", action="store_true", help="Verbose output",
    )
    args = parser.parse_args()

    if args.verbose:
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
        n_workers=12, processes=True, threads_per_worker=1, local_directory="/dev/shm"
    )
    client = Client(cluster)
    log.debug(client)

    main(
        polcat=args.polcat,
        client=client,
        data_dir=args.data_dir,
        do_update_cubes=args.update_cubes,
        do_convert_spectra=args.convert_spectra,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    cli()
