#!/usr/bin/env python3
"""Prepare files for CASDA upload"""
from email.policy import default
from astropy.io import fits
import numpy as np
import logging as log
import traceback

def find_spectra():
    pass


def convert_spectra(spectrum: str) -> str:
    """Convert a ascii spectrum to FITS

    Args:
        spectrum (str): Name of ASCII spectrum file

    Returns:
        str: Name of FITS spectrum file
    """    

    data = np.loadtxt(spectrum)
    # freq_Hz, I, Q, U, dI, dQ, dU
    try:
        log.info("> Trying [freq_Hz, I, Q, U, dI, dQ, dU]")
        (freqArr_Hz, IArr, QArr, UArr, dIArr, dQArr, dUArr) = data
        log.info("... success.")
    except Exception as e:
        log.info("...failed.")
        # freq_Hz, q, u, dq, du
        try:
            log.info("> Trying [freq_Hz, q, u,  dq, du]", end=' ')
            (freqArr_Hz, QArr, UArr, dQArr, dUArr) = data
            log.info("... success.")
            noStokesI = True
        except Exception as e:
            log.info("...failed.")
            log.debug(traceback.format_exc())
            raise e
    log.info("Successfully read in the Stokes spectra.")

    # If no Stokes I present, create a dummy spectrum = unity
    if noStokesI:
        log.warning("No Stokes I data in use.")
        IArr = np.ones_like(QArr)
        dIArr = np.zeros_like(QArr)

    # Create the FITS file
    hdu = fits.PrimaryHDU()
    hdu.header['NAXIS'] = 3
    hdu.header['NAXIS1'] = len(freqArr_Hz)
    hdu.header['NAXIS2'] = 2
    hdu.header['NAXIS3'] = 1
    hdu.header['CTYPE1'] = 'FREQ'
    hdu.header['CTYPE2'] = 'STOKES'
    hdu.header['CTYPE3'] = 'STOKES'
    hdu.header['CRVAL1'] = freqArr_Hz[0]
    hdu.header['CRVAL2'] = 1
    hdu.header['CRVAL3'] = 1
    hdu.header['CRPIX1'] = 1
    hdu.header['CRPIX2'] = 1
    hdu.header['CRPIX3'] = 1
    hdu.header['CDELT1'] = freqArr_Hz[1] - freqArr_Hz[0]
    hdu.header['CDELT2'] = 1
    hdu.header['CDELT3'] = 1
    hdu.header['CUNIT1'] = 'Hz'
    hdu.header['CUNIT2'] = 'STOKES'
    hdu.header['CUNIT3'] = 'STOKES'
    hdu.header['BUNIT'] = 'Jy'
    hdu.header['OBJECT'] = 'Spectrum'


def update_cubes():
    pass


def find_cubes():
    pass


def main(do_update_cubes: bool = False, do_covert_spectra: bool = False):
    """Main function"""

    if do_update_cubes:
        cubes = find_cubes()
        for cube in cubes:
            update_cubes(cube)
    if do_covert_spectra:
        spectra = find_spectra()
        for spectrum in spectra:
            convert_spectra(spectrum)


def cli():
    import argparse
    parser = argparse.ArgumentParser(
        description="Prepare files for CASDA upload")
    parser.add_argument("--update-cubes", action="store_true",
                        help="Update cubes", default=False)
    parser.add_argument("--convert-spectra", action="store_true",
                        help="Convert spectra", default=False)
    args = parser.parse_args()
    main(
        do_update_cubes=args.update_cubes,
        do_covert_spectra=args.convert_spectra
    )


if __name__ == "__main__":
    cli()
