#!/usr/bin/env python
from spiceracs.utils import gettable, tmatchn, tmatchtwo
import subprocess
import multiprocessing
import schwimmbad
import shlex
from tqdm import tqdm, trange
from glob import glob
import warnings
import os
import functools
print = functools.partial(print, flush=True)


def dobane(filename, n_cores=1, verbose=True):
    """Run BANE on moment map.

    Args:
        filename (str): FITS file to run BANE on.

    Kwargs:
        n_cores (int): Number of cores to run BANE with.
        verbose (bool): Print out messages.

    """
    #if verbose:
    #    print(f'Running BANE on {filename}')
    command = f"BANE {filename} --cores={n_cores}"
    command = shlex.split(command)
    #try:
    proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')


def doaegean(filename, moment, stoke, n_cores=1, verbose=True):
    """Run Aegean on moment map.

    Note: Assumes BANE has been run previously.

    Args:
        filename (str): FITS file to run BANE on.
        moment (str): Which 'Faraday' moment is contained within the
            file. Can be 'mu' or 'sigma'
        stoke (str): Which Stokes is containted within the file. Can be
            'p', 'q' or 'u'.

    Kwargs:
        n_cores (int): Number of cores to run BANE with.
        verbose (bool): Print out messages.

    """
    if verbose:
        print(f'Running Aegean on {filename}...')
    tablename = f"{filename.replace('.fits', '.xml')},\
    {filename.replace('.fits', '.reg')}".replace("   ", "").replace(
        " ", ""
    )
    # Defaults
    INNERCLIP = '5'
    OUTERCLIP = '4'

    if stoke == 'p':
        INNERCLIP = '6'
        OUTERCLIP = '5'

    command = ['aegean', filename, '--autoload', '--island',
               '--seedclip', INNERCLIP, '--floodclip', OUTERCLIP,
               f'--cores={n_cores}', '--table', tablename]

    if (stoke == 'q' or stoke == 'u') and moment == 'mu':
        negative = '--negative'
        command.append(negative)

    proc = subprocess.run(command,
                          capture_output=(not verbose),
                          encoding="utf-8", check=True)


def getmoments(momdir='.', zero=False, verbose=True):
    """Get moment files.

    Kwargs:
        momdir (str): Directory containing moment maps.
        verbose (bool): Print out messages.

    Returns:
        moments (dict): Filenames of momemt maps.

    """

    moments = {}
    for stoke in ['p', 'q', 'u']:
        mu = glob(f'{momdir}/*.{stoke}.mu.fits')
        sigma = glob(f'{momdir}/*.{stoke}.sigma.fits')
        moments.update({
            f"{stoke}_mu": mu,
            f"{stoke}_sigma": sigma
        })

    if zero:
        zero_file = glob(f'{momdir}/*.mom0.fits')
        moments.update({
            "p_mom0": zero_file[0],
        })
    return moments


def baneloop(moments, n_cores=1, verbose=True):
    """Loop BANE over files.

    Args:
        moments (dict): Filenames of momemt maps.

    Kwargs:
        n_cores (int): Number of cores to run BANE with.
        verbose (bool): Print out messages.

    """
    for key in tqdm(
            moments,
            disable=(not verbose),
            desc='Running BANE'
    ):
        for filename in tqdm(moments[key]):
            dobane(filename, n_cores=n_cores, verbose=verbose)


def aegeanloop(moments, n_cores=1, verbose=True):
    """Loop Aegean over files.

    Args:
        moments (dict): Filenames of momemt maps.

    Kwargs:
        n_cores (int): Number of cores to run BANE with.
        verbose (bool): Print out messages.

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


def squishtables(catdir='.', component=False, verbose=True):
    """Merge tables with STILTS.
    """
    if catdir[-1] == '/':
        catdir = catdir[:-1]

    cattype = 'isle'
    if component is True:
        cattype = 'comp'

    catalogues = {}
    for stoke in ['p', 'q', 'u']:
        for moment in ['mu', 'sigma']:
            tab, name = gettable(
                catdir,
                f'{stoke}*{moment}*linmos_{cattype}',
                verbose=verbose
            )
            catalogues.update({
                f"{stoke}_{moment}": name
            })
    catnames = [catalogues[key] for key in catalogues.keys()]

    # Run STILTS -- xmatch pol tables
    inN = [catnames[0], catnames[1]]
    valuesN = ['ra dec', 'ra dec']
    tmatchtwo(inN, valuesN, join='1or2', out='temp1.xml', verbose=verbose)

    for i in trange(1, 4):
        inN = [f'temp{i}.xml', catnames[i+2]]
        tmatchtwo(inN, valuesN, out=f'temp{i+1}.xml')

    # Run STILTS -- xmatch with I cat
    #inN = ['temp6.xml', icat]
    #valuesN = ['ra dec', 'ra_deg_cont dec_deg_cont']
    #tmatchtwo(inN, out='final.xml', join='1and2')

    # clean up
    # for i in range(5):
    #    os.remove(f'temp{i}.xml')


def main(args, verbose=True):
    """Main script.
    """
    # Sort out args
    momdir = args.momdir
    if momdir[-1] == '/':
            momdir = momdir[:-1]
    tabledir = args.tabledir

    n_cores = args.n_cores
    if args.n_cores is None:
        n_cores = multiprocessing.cpu_count()

    # Read in data
    if verbose:
        print('Finding data...')
    moments = getmoments(momdir, zero=args.do_zero, verbose=verbose)

    if args.do_zero:
        dobane(moments["p_mom0"], n_cores=n_cores, verbose=verbose)
        doaegean(moments["p_mom0"], moment='mom0', stoke='p',
                 n_cores=n_cores, verbose=verbose)

    # Run BANE
    if args.do_bane:
        baneloop(moments, n_cores=n_cores, verbose=verbose)

    # Run Aegean
    if args.do_aegean:
        aegeanloop(moments, n_cores=n_cores, verbose=verbose)

    # Merge tables
    if args.do_merge:
        squishtables(momdir, verbose=verbose)

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

    parser.add_argument(
        "--bane",
        dest="do_bane",
        action="store_true",
        help="Run BANE [False]."
    )

    parser.add_argument(
        "--aegean",
        dest="do_aegean",
        action="store_true",
        help="Run Aegean [False]."
    )

    parser.add_argument(
        "--zero",
        dest="do_zero",
        action="store_true",
        help="Run on zeroth moment [False]."
    )

    parser.add_argument(
        "--merge",
        dest="do_merge",
        action="store_true",
        help="Merge output tables [False]."
    )

    # TO-DO
    #parser.add_argument("--loners", dest="loners", action="store_true",
    #                    help="Run on single component sources [False].")

    args = parser.parse_args()
    verbose = args.verbose

    main(args, verbose=verbose)


if __name__ == "__main__":
    cli()
