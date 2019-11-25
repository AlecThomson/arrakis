#!/usr/bin/env python
from spiceracs.utils import gettable
import subprocess
import multiprocessing
from tqdm import tqdm, trange
from glob import glob
import warnings

def dobane(filename, n_cores=1, verbose=True):
    """Run BANE on moment map.
    """
    if verbose:
        print(f'Running BANE on {filename}')
    command = f'BANE {filename} --cores={n_cores}'
    proc = subprocess.run(command, shell=True,
                          capture_output=(not verbose),
                          encoding="utf-8", check=True)


def doaegean(filename, moment, stoke, n_cores=1, verbose=True):
    """Run Aegean on moment map.
    """
    if verbose:
        print(f'Running Aegean on {filename}...')
    tablename = filename.replace('.fits', '.vot')
    # Defaults
    INNERCLIP = 5
    OUTERCLIP = 4
    negative = ''

    if stoke == 'p':
        INNERCLIP = 6
        OUTERCLIP = 5
    elif (stoke == 'q' or stoke == 'u') and moment == 'mu':
        negative = '--negative'

    command = f'aegean {filename} --autoload --seedclip {INNERCLIP} \
        --floodclip {OUTERCLIP} --cores={n_cores} \
            {negative} --table {tablename}'
    proc = subprocess.run(command, shell=True,
                          capture_output=(not verbose),
                          encoding="utf-8", check=True)


def getmoments(momdir, verbose=True):
    """Get moment files
    """
    for stoke in ['p', 'q', 'u']:
        mu = glob(f'{momdir}/*.{stoke}.*.mu.*.fits')
        sigma = glob(f'{momdir}/*.{stoke}.*.sigma.*.fits')

        moments = {
            f"{stoke}_mu": mu,
            f"{stoke}_sigma": sigma
        }
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
    for stoke in tqdm(
            ['p', 'q', 'u'],
            disable=(not verbose),
            desc='Running Aegean'
            ):
        for moment in ['mu', 'sigma']:
            doaegean(moments[f'{stoke}_{moment}'], stoke, moment,
                     n_cores=n_cores, verbose=verbose)

def squishtables(Verbose=True):
    command = f'java -jar spiceracs/stilts.jar tmatch2\
         ifmt1=votable ifmt2=votable \
             in1=p.mu_comp.vot in2=u.sigma_comp.vot \
                 join=1or2 matcher=sky params=10 \
                     values1="RA Dec" values2="RA Dec" \
                         omode=out out=comb.vot'
    
    proc = subprocess.run(command, shell=True,
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
