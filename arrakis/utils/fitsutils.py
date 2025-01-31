#!/usr/bin/env python
"""FITS utilities"""

from __future__ import annotations

import warnings
from glob import glob
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS
from spectral_cube import SpectralCube
from spectral_cube.utils import SpectralCubeWarning

from arrakis.logger import logger
from arrakis.utils.io import gettable
from FRion.correct import find_freq_axis

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)


def head2dict(h: fits.Header) -> dict[str, Any]:
    """Convert FITS header to a dict.

    Writes a cutout, as stored in source_dict, to disk. The file location
    should already be specified in source_dict. This format is intended
    for parallel use with pool.map syntax.

    Args:
        h: An astropy FITS header.

    Returns:
        data (dict): The FITS head converted to a dict.

    """
    data = {}
    for c in h.__dict__["_cards"]:
        if c[0] == "":
            continue
        data[c[0]] = c[1]
    return data


def fix_header(cutout_header: fits.Header, original_header: fits.Header) -> fits.Header:
    """Make cutout header the same as original header

    Args:
        cutout_header (fits.Header): Cutout header
        original_header (fits.Header): Original header

    Returns:
        fits.Header: Fixed header
    """
    axis_cut = find_freq_axis(cutout_header)
    axis_orig = find_freq_axis(original_header)
    fixed_header = cutout_header.copy()
    if axis_cut != axis_orig:
        for key, val in cutout_header.items():
            if key[-1] == str(axis_cut):
                fixed_header[f"{key[:-1]}{axis_orig}"] = val
                fixed_header[key] = original_header[key]

    return fixed_header


def getfreq(
    cube: str | Path,
    outdir: Path | None = None,
    filename: str | Path | None = None,
) -> u.Quantity | tuple[u.Quantity, Path]:
    """Get list of frequencies from FITS data.

    Gets the frequency list from a given cube. Can optionally save
    frequency list to disk.

    Args:
        cube (str): File to get spectral axis from.

    Kwargs:
        outdir (str): Where to save the output file. If not given, data
            will not be saved to disk.

        filename (str): Name of frequency list file. Requires 'outdir'
            to also be specified.

        verbose (bool): Whether to print messages.

    Returns:
        freq (list): Frequencies of each channel in the input cube.

    """
    with fits.open(cube, memmap=True, mode="denywrite") as hdulist:
        hdu = hdulist[0]
        hdr = hdu.header
        data = hdu.data

    # Two problems. The default 'UTC' stored in 'TIMESYS' is
    # incompatible with the TIME_SCALE checks in astropy.
    # Deleting or coverting to lower case fixes it. Second
    # problem, the OBSGEO keywords prompts astropy to apply
    # a velocity correction, but no SPECSYS has been defined.
    for k in ["TIMESYS", "OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z"]:
        if k in hdr:
            del hdr[k]

    wcs = WCS(hdr)
    freq: u.Quantity = wcs.spectral.pixel_to_world(np.arange(data.shape[0]))

    # Write to file if outdir is specified
    if outdir is None:
        return freq

    outfile = outdir / filename if filename is not None else outdir / "frequencies.txt"
    logger.info(f"Saving to {outfile}")
    np.savetxt(outfile, np.array(freq))
    return freq, outfile


def getdata(cubedir="./", tabledir="./", mapdata=None, verbose=True):
    """Get the spectral and source-finding data.

    Args:
        cubedir: Directory containing data cubes in FITS format.
        tabledir: Directory containing Selavy results.
        mapdata: 2D FITS image which corresponds to Selavy table.

    Kwargs:
        verbose (bool): Whether to print messages.

    Returns:
        datadict (dict): Dictionary of necessary astropy tables and
            Spectral cubes.

    """
    if cubedir[-1] == "/":
        cubedir = cubedir[:-1]

    if tabledir[-1] == "/":
        tabledir = tabledir[:-1]
    # Glob out the necessary files
    # Data cubes
    icubes = glob(f"{cubedir}/image.restored.i.*contcube*linmos.fits")
    qcubes = glob(f"{cubedir}/image.restored.q.*contcube*linmos.fits")
    ucubes = glob(f"{cubedir}/image.restored.u.*contcube*linmos.fits")
    vcubes = glob(f"{cubedir}/image.restored.v.*contcube*linmos.fits")

    cubes = [icubes, qcubes, ucubes, vcubes]
    # Selavy images
    selavyfits = mapdata
    # Get selvay data from VOTab
    i_tab, voisle = gettable(tabledir, "islands", verbose=verbose)  # Selvay VOTab
    components, tablename = gettable(tabledir, "components", verbose=verbose)

    logger.info(f"Getting spectral data from: {cubes}\n")
    logger.info(f"Getting source location data from: {selavyfits}\n")

    # Read data using Spectral cube
    i_taylor = SpectralCube.read(selavyfits, mode="denywrite")
    wcs_taylor = WCS(i_taylor.header)
    i_cube = SpectralCube.read(icubes[0], mode="denywrite")
    wcs_cube = WCS(i_cube.header)
    q_cube = SpectralCube.read(qcubes[0], mode="denywrite")
    u_cube = SpectralCube.read(ucubes[0], mode="denywrite")
    if len(vcubes) != 0:
        v_cube = SpectralCube.read(vcubes[0], mode="denywrite")
    else:
        v_cube = None
    # Mask out using Stokes I == 0 -- seems to be the current fill value
    mask = ~(i_cube == 0 * u.jansky / u.beam)
    i_cube = i_cube.with_mask(mask)
    mask = ~(q_cube == 0 * u.jansky / u.beam)
    q_cube = q_cube.with_mask(mask)
    mask = ~(u_cube == 0 * u.jansky / u.beam)
    u_cube = u_cube.with_mask(mask)

    datadict = {
        "i_tab": i_tab,
        "i_tab_comp": components,
        "i_taylor": i_taylor,
        "wcs_taylor": wcs_taylor,
        "wcs_cube": wcs_cube,
        "i_cube": i_cube,
        "q_cube": q_cube,
        "u_cube": u_cube,
        "v_cube": v_cube,
        "i_file": icubes[0],
        "q_file": qcubes[0],
        "u_file": ucubes[0],
        "v_file": vcubes[0],
    }

    return datadict
