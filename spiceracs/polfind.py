#!/usr/bin/env python
from spiceracs.utils import gettable
import subprocess
import multiprocessing
from tqdm import tqdm, trange
from glob import glob
import warnings
import functools
print = functools.partial(print, flush=True)


def dobane(filename, n_cores=1, verbose=True):
    """Run BANE on moment map.
    """
    if verbose:
        print(f'Running BANE on {filename}')
    command = ['BANE', filename, f'--cores={n_cores}']
    proc = subprocess.run(command,
                          capture_output=(not verbose),
                          encoding="utf-8", check=True)


def doaegean(filename, moment, stoke, n_cores=1, verbose=True):
    """Run Aegean on moment map.
    """
    if verbose:
        print(f'Running Aegean on {filename}...')
    tablename = f"{filename.replace('.fits', '.vot')},\
    {filename.replace('.fits', '.reg')}".replace("   ", "").replace(
        " ", ""
    )
    # Defaults
    INNERCLIP = 5
    OUTERCLIP = 4
    negative = ''

    if stoke == 'p':
        INNERCLIP = 6
        OUTERCLIP = 5
    elif (stoke == 'q' or stoke == 'u') and moment == 'mu':
        negative = '--negative'

    command = ['aegean', filename, '--autoload', '--island',
               f'--seedclip {INNERCLIP}', f'--floodclip {OUTERCLIP}',
               f'--cores={n_cores}', negative, f'--table {tablename}']
    proc = subprocess.run(command,
                          capture_output=(not verbose),
                          encoding="utf-8", check=True)


def getmoments(momdir, verbose=True):
    """Get moment files
    """
    moments = {}
    for stoke in ['p', 'q', 'u']:
        mu = glob(f'{momdir}/*.{stoke}.*.mu.*linmos.fits')
        sigma = glob(f'{momdir}/*.{stoke}.*.sigma.*.linmos.fits')
        moments.update({
            f"{stoke}_mu": mu[0],
            f"{stoke}_sigma": sigma[0]
        })
    return moments


def baneloop(moments, n_cores=1, verbose=True):
    """Loop BANE over files.
    """
    for key in tqdm(
            moments,
            disable=(not verbose),
            desc='Running BANE'
    ):
        dobane(moments[key], n_cores=n_cores, verbose=verbose)


def aegeanloop(moments, n_cores=1, verbose=True):
    """Loop Aegean over files.
    """
    with tqdm(
        total=6,
        disable=(not verbose),
        desc='Running Aegean'
    ) as pbar:
        for stoke in ['p', 'q', 'u']:
            for moment in ['mu', 'sigma']:
                doaegean(moments[f'{stoke}_{moment}'], moment, stoke,
                         n_cores=n_cores, verbose=verbose)
                pbar.update(1)


def squishtables(verbose=True):
    stiltspath = Path(os.path.realpath(__file__)
                      ).parent.parent/"stilts"/"stilts.jar"
    command = ['java','-jar', stiltspath,'-h']
    proc = subprocess.run(command,
                        capture_output=(not verbose),
                        encoding="utf-8", check=True)


def main(args, verbose=True):
    """Main script.
    """
    # Sort out args
    momdir = args.momdir
    tabledir = args.tabledir
    if momdir[-1] == '/':
        momdir = momdir[:-1]

    if tabledir[-1] == '/':
        tabledir = tabledir[:-1]

    n_cores = args.n_cores
    if args.n_cores is None:
        n_cores = multiprocessing.cpu_count()

    # Read in data
    if verbose:
        print('Finding data...')
    moments = getmoments(momdir, verbose=verbose)

    # Run BANE
    baneloop(moments, n_cores=n_cores, verbose=verbose)

    # Run Aegean
    aegeanloop(moments, n_cores=n_cores, verbose=verbose)

    # Merge tables
    squishtables()

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
    SPICE-RACS Stage 4:
    Polarized source-finding.


    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        'momdir',
        metavar='momdir',
        type=str,
        help='Directory containing moment maps.')

    parser.add_argument(
        'tabledir',
        metavar='tabledir',
        type=str,
        help='Directory containing Selavy results.')

    parser.add_argument(
        "-v",
        dest="verbose",
        action="store_true",
        help="Verbose output [False]."
    )

    parser.add_argument(
        "--ncores",
        dest="n_cores",
        type=int, help="Number of processes [use all available].")

    args = parser.parse_args()
    verbose = args.verbose

    main(args, verbose=verbose)


if __name__ == "__main__":
    cli()
