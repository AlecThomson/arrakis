#!/usr/bin/env python
from spiceracs.cutout import getdata
import warnings
import numpy as np
from astropy.io import fits
from tqdm import tqdm, trange

class moments:
    def __init__(self, cube, pool, verbose=True):
        self.cube = cube
        self.verbose = verbose
        self.pool = pool

    def mu(self, outdir='.', stoke=''):
        """
        Mean moment - freq axis first
        """
        momfile = f'{outdir}/moment.mu.{stoke}.fits'
        blank = self.cube[0]
        blank.write(momfile, overwrite=True, format='fits')
        
        if self.verbose:
            print(f'Writing to {momfile}...')
            print(f'Looping over y-axis to save memory...')
        with fits.open(momfile, mode='update', memmap=True) as outfh:
            for i in trange(
                    self.cube.shape[1],
                    desc='Looping over y-axis'
                ):
                yav = self.cube[:,i,:].mean(axis=0)
                outfh[0].data[i,:] = yav
                outfh.flush()
        with fits.open(momfile, mode='denywrite') as outfh:
            print(outfh[0].data)

    def muworker(args):
        i, strip, outfh = args
        yav = strip.mean(axis=0)
        outfh[0].data[i,:] = yav
        outfh.flush()
        
    def sigma(self):
        """
        Standard deviation moment - freq axis first
        """
        n = self.cube.shape[0]
        sigma = np.sqrt((1/(n-1)) * np.nansum((self.cube - mu)**2, axis=0))
        return sigma


def momentloop(datadict, outdir='.', verbose=True):
    for stokes in ['q', 'u']:
        if verbose:
            print(f'Making moments for Stokes {stokes}...')
        mom = moments(datadict[f'{stokes}_cube'], verbose=verbose)
        outfile = f'{outdir}/moment.{stokes}.fits'
        mom.mu(outdir=outdir, stoke=stokes)
        


def main(args, verbose=True):
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
    
    # Compute moments
    momentloop(datadict, outdir, verbose=verbose)

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
    Polarized source-finding.

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

    parser.add_argument(
        "-d",
        dest="dryrun",
        action="store_true",
        help="Do a dry-run [False]."
    )

    parser.add_argument(
        "-s",
        dest="getsize",
        action="store_true",
        help="Estimate size of cutouts [False]."
    )

    parser.add_argument(
        "-m",
        dest="database",
        action="store_true",
        help="Add data to MongoDB [False]."
    )

    args = parser.parse_args()
    verbose = args.verbose

    if args.database:
        if verbose:
            print('Testing MongoDB connection...')
        client = pymongo.MongoClient()  # default connection (ie, local)
        try:
            client.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")
        else:
            if verbose:
                print('MongoDB connection succesful!')
        client.close()

    main(args, verbose=verbose)

if __name__ == "__main__":
    cli()
