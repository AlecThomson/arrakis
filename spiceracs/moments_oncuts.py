#!/usr/bin/env python
from spiceracs.utils import getfreq
import numpy as np
import os
import pymongo
import sys
import subprocess
import time
from tqdm import tqdm, trange
import warnings
from RMtools_1D import do_RMsynth_1D
from RMtools_3D import do_RMsynth_3D
from spectral_cube import SpectralCube
from astropy.io import fits
import astropy.units as u
import functools
import psutil
print = functools.partial(print,f'[{psutil.Process().cpu_num()}]', flush=True)


def moment_worker(args):
    """Make moments of cutouts
    """
    i, clargs, outdir, verbose = args

    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

    if clargs.pol and not clargs.unres:
        myquery = {"polarized": True}
    elif clargs.unres and not clargs.pol:
        myquery = {"resolved": False}
    elif clargs.pol and clargs.unres:
        myquery = {"$and": [{"resolved": False}, {"polarized": True}]}

    elif clargs.pol and not clargs.loners:
        myquery = {"polarized": True}
    elif clargs.loners and not clargs.pol:
        myquery = {"n_components": 1}
    elif clargs.pol and clargs.loners:
        myquery = {"$and": [{"n_components": 1}, {"polarized": True}]}
    else:
        myquery = {}

    doc = mycol.find(myquery).sort("flux_peak", -1)
    iname = doc[i]['island_name']

    for stokes in ['p', 'q', 'u']:
        infile = doc[i][f'{stokes}_file']
        mufile = infile.replace('.fits', '.mu.fits')
        sigfile = infile.replace('.fits', '.sigma.fits')

        data = SpectralCube.read(f"{outdir}/{infile}")

        mu = data.mean(axis=0)
        sigma = data.std(axis=0, ddof=1)

        mu.write(f"{outdir}/{mufile}", format='fits', overwrite=True)
        sigma.write(f"{outdir}/{sigfile}", format='fits', overwrite=True)

        if clargs.database:
            myquery = {"island_name": iname}
            newvalues = {"$set": {f"{stokes}_mu_file": mufile}}
            mycol.update_one(myquery, newvalues)
            newvalues = {"$set": {f"{stokes}_sigma_file": sigfile}}
            mycol.update_one(myquery, newvalues)


def moment_worker_i(args):
    """Make moments of cutouts
    """
    i, clargs, outdir, verbose = args

    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

    if clargs.pol and not clargs.unres:
        myquery = {"polarized": True}
    elif clargs.unres and not clargs.pol:
        myquery = {"resolved": False}
    elif clargs.pol and clargs.unres:
        myquery = {"$and": [{"resolved": False}, {"polarized": True}]}

    elif clargs.pol and not clargs.loners:
        myquery = {"polarized": True}
    elif clargs.loners and not clargs.pol:
        myquery = {"n_components": 1}
    elif clargs.pol and clargs.loners:
        myquery = {"$and": [{"n_components": 1}, {"polarized": True}]}
    else:
        myquery = {}

    doc = mycol.find(myquery).sort("flux_peak", -1)
    iname = doc[i]['island_name']

    stokes = 'i'
    infile = doc[i][f'{stokes}_file']
    mufile = infile.replace('.fits', '.mu.fits')

    data = SpectralCube.read(f"{outdir}/{infile}")

    mu = data.mean(axis=0)

    mu.write(f"{outdir}/{mufile}", format='fits', overwrite=True)

    if clargs.database:
        myquery = {"island_name": iname}
        newvalues = {"$set": {f"{stokes}_mu_file": mufile}}
        mycol.update_one(myquery, newvalues)

def makepiworker(args):
    """Make PI of cutouts
    """
    i, clargs, outdir, verbose = args

    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

    if clargs.pol and not clargs.unres:
        myquery = {"polarized": True}
    elif clargs.unres and not clargs.pol:
        myquery = {"resolved": False}
    elif clargs.pol and clargs.unres:
        myquery = {"$and": [{"resolved": False}, {"polarized": True}]}

    elif clargs.pol and not clargs.loners:
        myquery = {"polarized": True}
    elif clargs.loners and not clargs.pol:
        myquery = {"n_components": 1}
    elif clargs.pol and clargs.loners:
        myquery = {"$and": [{"n_components": 1}, {"polarized": True}]}
    else:
        myquery = {}

    doc = mycol.find(myquery).sort("flux_peak", -1)
    iname = doc[i]['island_name']

    qfile = doc[i][f'q_file']
    ufile = doc[i][f'u_file']
    outfile = qfile.replace('.q.', '.p.')

    qdata = SpectralCube.read(f"{outdir}/{qfile}")
    udata = SpectralCube.read(f"{outdir}/{ufile}")

    pdata = (qdata**2 + udata**2)**0.5

    pdata.write(f'{outdir}/{outfile}', format='fits', overwrite=True)

    if clargs.database:
        myquery = {"island_name": iname}
        newvalues = {"$set": {"p_file": outfile}}
        mycol.update_one(myquery, newvalues)


def zero_worker(args):
    """Make Zeroth Faraday moment
    """
    i, clargs, outdir, verbose = args

    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

    doc = mycol.find().sort("flux_peak", -1)
    iname = doc[i]['island_name']

    qfile = doc[i][f'q_file']
    ufile = doc[i][f'u_file']
    outfile = qfile.replace('.q.', '.p.').replace('.fits', '.mom0.fits')

    qdata = SpectralCube.read(f"{outdir}/{qfile}")
    udata = SpectralCube.read(f"{outdir}/{ufile}")

    freq = np.array(getfreq(qdata))

    dataArr = do_RMsynth_3D.run_rmsynth(
                np.array(qdata),
                np.array(udata),
                freq, phiMax_radm2=clargs.phiMax_radm2, nSamples=5,
                weightType="uniform", fitRMSF=False, nBits=32,
                verbose=False, not_rmsf=True
            )
    FDFcube, phiArr_radm2, lam0Sq_m2, lambdaSqArr_m2 = dataArr
    dphi = np.diff(phiArr_radm2)[0]
    mom0 = np.nansum(abs(FDFcube)*dphi, axis=0)
    newf = fits.PrimaryHDU()
    newf.data = mom0
    newf.header = qdata.header
    newf.header.update(qdata[0].wcs.to_header())
    newf.writeto(f'{outdir}/{outfile}', overwrite=True)


    if clargs.database:
        myquery = {"island_name": iname}
        newvalues = {"$set": {"p_mom0_file": outfile}}
        mycol.update_one(myquery, newvalues)

def main(pool, args, verbose=False):
    """Main script
    """
    outdir = args.outdir
    if outdir[-1] == '/':
        outdir = outdir[:-1]
    outdir = f'{outdir}/cutouts'
    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

        # Basic querey
    if args.pol and not args.unres:
        myquery = {"polarized": True}
    elif args.unres and not args.pol:
        myquery = {"resolved": False}
    elif args.pol and args.unres:
        myquery = {"$and": [{"resolved": False}, {"polarized": True}]}

    elif args.pol and not args.loners:
        myquery = {"polarized": True}
    elif args.loners and not args.pol:
        myquery = {"n_components": 1}
    elif args.pol and args.loners:
        myquery = {"$and": [{"n_components": 1}, {"polarized": True}]}

    else:
        myquery = {}

    mydoc = mycol.find(myquery).sort("flux_peak", -1)
    count = mycol.count_documents({})

    if args.limit is not None:
        count = args.limit

    if verbose:
        print(f'Making moments of {count} sources')

    inputs = [[i, args, outdir, verbose] for i in range(count)]
    if (pool.__class__.__name__ is 'MPIPool' or
            pool.__class__.__name__ is 'SerialPool'):
        if args.picube:
            # Make PI
            if verbose:
                print('Making PI...')
            tic = time.perf_counter()
            list(pool.map(makepiworker, inputs))
            toc = time.perf_counter()
            if verbose:
                print(f'Time taken was {toc - tic}s')
        if args.farnes:
            # Make moments
            if verbose:
                print('Making moments...')
            tic = time.perf_counter()
            list(pool.map(moment_worker, inputs))
            toc = time.perf_counter()
            if verbose:
                print(f'Time taken was {toc - tic}s')
        if args.zero:
            # Make moments
            if verbose:
                print('Making 0th moment...')
            tic = time.perf_counter()
            list(pool.map(zero_worker, inputs))
            toc = time.perf_counter()
            if verbose:
                print(f'Time taken was {toc - tic}s')

        if args.on_i:
            # Make moments
            if verbose:
                print('Making Stokes I moment...')
            tic = time.perf_counter()
            list(pool.map(moment_worker_i, inputs))
            toc = time.perf_counter()
            if verbose:
                print(f'Time taken was {toc - tic}s')
        

    elif pool.__class__.__name__ is 'MultiPool':
        if args.picube:
            list(tqdm(
                pool.imap_unordered(makepiworker, inputs),
                total=count,
                desc='Making PI',
                disable=(not verbose)
            )
            )
        if args.farnes:
            list(tqdm(
                pool.imap_unordered(moment_worker, inputs),
                total=count,
                desc='Making moments',
                disable=(not verbose)
            )
            )
        if args.zero:
            list(tqdm(
                pool.imap_unordered(zero_worker, inputs),
                total=count,
                desc='Making 0th moment',
                disable=(not verbose)
            )
            )

        if args.on_i:
            list(tqdm(
                pool.imap_unordered(moment_worker_i, inputs),
                total=count,
                desc='Making Stokes I moments',
                disable=(not verbose)
            )
            )

    pool.close()

    if verbose:
        print('Done!')


def cli():
    """Command-line interface
    """
    import argparse
    import schwimmbad
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
    from astropy.io.fits.verify import VerifyWarning
    warnings.simplefilter('ignore', category=VerifyWarning)
    warnings.filterwarnings(
        "ignore", message="Degrees of freedom <= 0 for slice.")
    warnings.filterwarnings(
        "ignore", message="Mean of empty slice")

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

    # Help string to be shown using the -h option
    descStr = f"""
    {logostr}
    SPICE-RACS Stage 3:
    Make moments from cubelets.

    Note: Runs on brightest sources first.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        'outdir',
        metavar='outdir',
        type=str,
        help='Directory containing cutouts (in subdir outdir/cutouts).')

    parser.add_argument("-v", dest="verbose", action="store_true",
                        help="verbose output [False].")

    parser.add_argument(
        "-m",
        dest="database",
        action="store_true",
        help="Add data to MongoDB [False]."
    )

    parser.add_argument("--limit", dest="limit", default=None,
                        type=int, help="Limit number of sources [All].")

    parser.add_argument(
        "--farnes",
        dest="farnes",
        action="store_true",
        help="Make Farnes (2018) moments [False]."
    )

    parser.add_argument(
        "--zero",
        dest="zero",
        action="store_true",
        help="Make Zeroth moment using RM synthesis [False].")

    parser.add_argument("-l", dest="phiMax_radm2", type=float, default=None,
                        help="Absolute max Faraday depth sampled (for mom0) [Auto].")

    parser.add_argument(
        "--picube",
        dest="picube",
        action="store_true",
        help="Make PI cubes [False].")

    parser.add_argument(
        "-i", "--stokesi",
        dest="on_i",
        action="store_true",
        help="Make Stokes I moments [False].")

    parser.add_argument("--pol", dest="pol", action="store_true",
                        help="Run on polarized sources [False].")

    parser.add_argument("--unres", dest="unres", action="store_true",
                        help="Run on unresolved sources [False].")

    parser.add_argument("--loners", dest="loners", action="store_true",
                        help="Run on single component sources [False].")

    group = parser.add_mutually_exclusive_group()

    group.add_argument("--ncores", dest="n_cores", default=1,
                       type=int, help="Number of processes (uses multiprocessing).")
    group.add_argument("--mpi", dest="mpi", default=False,
                       action="store_true", help="Run with MPI.")

    args = parser.parse_args()
    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)

    verbose = args.verbose

    if args.mpi:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)

    if verbose:
        print(f"Using pool: {pool.__class__.__name__}")

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

    main(pool, args, verbose=verbose)


if __name__ == "__main__":
    cli()
