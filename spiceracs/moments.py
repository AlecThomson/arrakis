#!/usr/bin/env python
from spiceracs.utils import getdata, copyfile
from spectral_cube import SpectralCube
import warnings
import numpy as np
from astropy.io import fits
from tqdm import tqdm, trange
import os
import pdb
import pymongo


class moments:
    """Methods for produce moment maps.
    """

    def __init__(self, cube, verbose=True):
        self.cube = cube
        self.verbose = verbose

    def mu(self, outdir='.', stoke=''):
        """Mean moment - freq axis first
        """
        momfile = f'{outdir}/moment.mu.{stoke}.fits'
        blank = self.cube[0]
        blank.write(momfile, overwrite=True, format='fits')

        if self.verbose:
            print(f'Writing to {momfile}...')
        with fits.open(momfile, mode='update', memmap=True) as outfh:
            for i in trange(
                self.cube.shape[1],
                desc='Looping over y-axis to save memory',
                disable=(not self.verbose)
            ):
                yav = self.cube[:, i, :].mean(axis=0)
                outfh[0].data[i, :] = yav
                outfh.flush()

    def sigma(self, outdir='.', stoke=''):
        """
        Standard deviation moment - freq axis first
        """
        momfile = f'{outdir}/moment.sigma.{stoke}.fits'
        blank = self.cube[0]
        blank.write(momfile, overwrite=True, format='fits')

        if self.verbose:
            print(f'Writing to {momfile}...')
        with fits.open(momfile, mode='update', memmap=True) as outfh:
            for i in trange(
                self.cube.shape[1],
                desc='Looping over y-axis to save memory',
                disable=(not self.verbose)
            ):
                ystd = self.cube[:, i, :].std(axis=0, ddof=1)
                outfh[0].data[i, :] = ystd
                outfh.flush()


def momentloop(datadict, outdir='.', verbose=True):
    """Produce moment maps.

    Args:
        datadict (dict): Dictionary containing tables and cubes.

    Kwargs:
        outdir (str): Directory to save data to.
        verbose (bool): Print out messages.

    Returns:

    """
    for stokes in ['p', 'q', 'u']:
        if verbose:
            print(f'Making moments for Stokes {stokes}...')
        mom = moments(datadict[f'{stokes}_cube'], verbose=verbose)
        outfile = f'{outdir}/moment.{stokes}.fits'
        if verbose:
            print(f'Making mean...')
        mom.mu(outdir=outdir, stoke=stokes)
        if verbose:
            print(f'Making std...')
        mom.sigma(outdir=outdir, stoke=stokes)


def makepi(datadict, verbose=True):
    """Make a polarized itensity cube.

    Note: PI cube will be saved to same dir as other cubes.

    Args:
        datadict (dict): Dictionary containing tables and cubes.

    Kwargs:
        verbose (bool): Print out messages.

    Returns:
        datadict (dict): Dictionary containing tables and cubes, now
        updated with PI data.

    """
    pifile = datadict['i_file'].replace('.i.', '.p.')
    if verbose:
        print(f'Writing to {pifile}...')
    try:
        with open(pifile, 'rb') as fh:
            pass
    except FileNotFoundError:
        copyfile(datadict['i_file'], pifile, verbose=verbose)

    with fits.open(pifile, mode='update', memmap=True) as outfh:
        for i in trange(
            datadict['i_cube'].shape[1],
            desc='Looping over y-axis to save memory',
            disable=(not verbose)
        ):
            ypi = np.hypot(
                np.expand_dims(datadict['q_cube'][:, i, :], axis=1),
                np.expand_dims(datadict['u_cube'][:, i, :], axis=1)
            )
            outfh[0].data[:, :, i, :] = ypi
            outfh.flush()

    p_cube = SpectralCube.read(pifile, mode='denywrite')
    datadict['p_file'] = pifile
    datadict['p_cube'] = p_cube

    return datadict


def main(args, verbose=True):
    """Main script.
    """
    # Sort out args
    cubedir = args.cubedir
    tabledir = args.tabledir
    if cubedir[-1] == '/':
        cubedir = cubedir[:-1]

    if tabledir[-1] == '/':
        tabledir = tabledir[:-1]

    outdir = args.outdir
    if outdir[-1] == '/':
        outdir = outdir[:-1]
    # Read in data
    if verbose:
        print('Reading data...')
    datadict = getdata(cubedir, tabledir, verbose=verbose)

    # Make PI cube
    if verbose:
        print('Making polarized intensity cube...')
    datadict = makepi(datadict, verbose=verbose)

    # Compute moments
    if verbose:
        print('Making moment maps...')
    momentloop(datadict, outdir=outdir, verbose=verbose)

    if verbose:
        print('Done!')


def cli():
    """Command-line interface.
    """
    import argparse
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
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

    descStr = f"""
    {logostr}
    SPICE-RACS Stage 3:
    Produce 'Faraday' moments.

    Note: A PI cube will also be produced, and saved with other cubes.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        'cubedir',
        metavar='cubedir',
        type=str,
        help='Directory containing data cubes in FITS format.')

    parser.add_argument(
        'tabledir',
        metavar='tabledir',
        type=str,
        help='Directory containing Selavy results.')

    parser.add_argument(
        'outdir',
        metavar='outdir',
        type=str,
        help='Directory to store moment maps.')

    parser.add_argument(
        "-v",
        dest="verbose",
        action="store_true",
        help="Verbose output [False]."
    )

    args = parser.parse_args()
    verbose = args.verbose

    main(args, verbose=verbose)


if __name__ == "__main__":
    cli()
