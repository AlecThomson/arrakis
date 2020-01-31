from spiceracs.utils import getfreq, MyEncoder
import json
import numpy as np
import os
import pymongo
import sys
import subprocess
import time
from tqdm import tqdm, trange
import warnings
from RMtools_1D import do_RMclean_1D
from RMtools_3D import do_RMclean_3D
from spectral_cube import SpectralCube
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from RMutils.util_misc import create_frac_spectra


def rmclean1d(args):
    i, clargs, outdir, verbose = args

    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

    # Basic querey
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
    prefix = f'{outdir}/{iname}'

    fdfFile = prefix + "_FDFdirty.dat"
    rmsfFile = prefix + "_RMSF.dat"
    weightFile = prefix + "_weight.dat"
    rmSynthFile = prefix + "_RMsynth.json"
    # Sanity checks
    for f in [weightFile, fdfFile, rmsfFile, rmSynthFile]:
        if not os.path.exists(f):
            print("File does not exist: '{:}'.".format(f), end=' ')
            sys.exit()
    nBits = 32
    mDictS, aDict = do_RMclean_1D.readFiles(
        fdfFile, rmsfFile, weightFile, rmSynthFile, nBits)
    # Run RM-CLEAN on the spectrum
    outdict, arrdict = do_RMclean_1D.run_rmclean(mDictS=mDictS,
                                                 aDict=aDict,
                                                 cutoff=clargs.cutoff,
                                                 maxIter=clargs.maxIter,
                                                 gain=clargs.gain,
                                                 nBits=nBits,
                                                 showPlots=clargs.showPlots,
                                                 verbose=clargs.rm_verbose)

    # Save output
    do_RMclean_1D.saveOutput(outdict,
                             arrdict,
                             prefixOut=prefix,
                             outDir=outdir,
                             verbose=clargs.rm_verbose)

    if clargs.database:
        # Load into Mongo
        myquery = {"island_name": iname}

        newvalues = {"$set": {"rmclean1d": True}}
        mycol.update_one(myquery, newvalues)

        newvalues = {"$set": {"rm_summary": mDictS}}
        mycol.update_one(myquery, newvalues)


def rmclean3d(args):
    i, clargs, outdir, verbose = args

    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

    # Basic querey
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
    qfile = f"{outdir}/{doc[i]['q_file']}"
    ufile = f"{outdir}/{doc[i]['u_file']}"
    vfile = f"{outdir}/{doc[i]['v_file']}"

    
    cleanFDF, ccArr, iterCountArr, residFDF, headtemp = do_RMclean_3D.run_rmclean(
        fitsFDF=f"{outdir}/{iname}FDF_real_dirty.fits",
        fitsRMSF=f"{outdir}/{iname}RMSF_tot.fits",
        cutoff=clargs.cutoff,
        maxIter=clargs.maxIter,
        gain=clargs.gain,
        chunksize=None,
        nBits=32,
        verbose=clargs.rm_verbose)

    # Write results to disk
    do_RMclean_3D.writefits(
        cleanFDF,
        ccArr,
        iterCountArr,
        residFDF,
        headtemp,
        prefixOut=iname,
        outDir=outdir,
        write_separate_FDF=True,
        verbose=clargs.rm_verbose)


def main(pool, args, verbose=False):
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
    count = mycol.count_documents(myquery)

    if args.limit is not None:
        count = args.limit

    if verbose:
        print(f'Running RM-CLEAN on {count} sources')

    if args.dimension == '1d':
        inputs = [[i, args, outdir, verbose]
                  for i in range(count)]
        if (pool.__class__.__name__ is 'MPIPool' or
                pool.__class__.__name__ is 'SerialPool'):
            if verbose:
                print('Running 1D RM synth...')
            tic = time.perf_counter()
            list(pool.map(rmclean1d, inputs))
            toc = time.perf_counter()
            if verbose:
                print(f'Time taken was {toc - tic}s')

        elif pool.__class__.__name__ is 'MultiPool':
            list(tqdm(
                pool.imap_unordered(rmclean1d, inputs),
                total=count,
                desc='Running 1D RM synth',
                disable=(not verbose)
            )
            )

    elif args.dimension == '3d':
        inputs = [[i, args, outdir, verbose]
                  for i in range(count)]
        if (pool.__class__.__name__ is 'MPIPool' or
                pool.__class__.__name__ is 'SerialPool'):
            if verbose:
                print('Running 3D RM-CLEAN...')
            tic = time.perf_counter()
            list(pool.map(rmclean3d, inputs))
            toc = time.perf_counter()
            if verbose:
                print(f'Time taken was {toc - tic}s')

        elif pool.__class__.__name__ is 'MultiPool':
            list(tqdm(
                pool.imap_unordered(rmclean3d, inputs),
                total=count,
                desc='Running 3D RM-CLEAN',
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
    SPICE-RACS Stage 6:
    Run RM-CLEAN on cubelets.

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

    parser.add_argument("--dimension", dest="dimension", default="1d",
                        help="How many dimensions for RMsynth [1d] or '3d'.")

    parser.add_argument("-v", dest="verbose", action="store_true",
                        help="verbose output [False].")

    parser.add_argument(
        "-m",
        dest="database",
        action="store_true",
        help="Add data to MongoDB [False]."
    )

    parser.add_argument("--pol", dest="pol", action="store_true",
                        help="Run on polarized sources [False].")

    parser.add_argument("--unres", dest="unres", action="store_true",
                        help="Run on unresolved sources [False].")

    parser.add_argument("--validate", dest="validate", action="store_true",
                        help="Run on Stokes I [False].")

    parser.add_argument("--limit", dest="limit", default=None,
                        type=int, help="Limit number of sources [All].")

    parser.add_argument("--loners", dest="loners", action="store_true",
                        help="Run on single component sources [False].")

    # RM-tools args
    parser.add_argument("-c", dest="cutoff", type=float, default=-3,
                        help="CLEAN cutoff (+ve = absolute, -ve = sigma) [-3].")
    parser.add_argument("-n", dest="maxIter", type=int, default=1000,
                        help="maximum number of CLEAN iterations [1000].")
    parser.add_argument("-g", dest="gain", type=float, default=0.1,
                        help="CLEAN loop gain [0.1].")
    parser.add_argument("-p", dest="showPlots", action="store_true",
                        help="show the plots [False].")
    parser.add_argument("-rmv", dest="rm_verbose", action="store_true",
                        help="Verbose RM-CLEAN [False].")

    group = parser.add_mutually_exclusive_group()

    group.add_argument("--ncores", dest="n_cores", default=1,
                       type=int, help="Number of processes (uses multiprocessing).")
    group.add_argument("--mpi", dest="mpi", default=False,
                       action="store_true", help="Run with MPI.")

    args = parser.parse_args()

    if args.validate:
        pool = schwimmbad.SerialPool
    else:
        pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)

    verbose = args.verbose

    if args.mpi:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)

    if verbose:
        print(f"Using pool: {pool.__class__.__name__}")

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

    main(pool, args, verbose=verbose)


if __name__ == "__main__":
    cli()
