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
import functools
import psutil
print = functools.partial(print, f'[{psutil.Process().cpu_num()}]', flush=True)


def rmclean1d(args):
    comp_id, clargs, outdir, host, field, verbose = args

    with pymongo.MongoClient(host=host) as client:
        mydb = client['spiceracs']  # Create/open database
        isl_col = mydb['islands']  # Create/open collection
        comp_col = mydb['components']  # Create/open collection
        beams_col = mydb['beams']  # Create/open collection

    # Basic querey
    myquery = {"Component_ID": comp_id}

    doc = comp_col.find_one(myquery)

    iname = doc['Source_ID']
    beam = beams_col.find_one({'Source_ID': iname})

    if clargs.rm_verbose:
        print(f'Working on {comp_id}')
    try:
        cname = comp_id
        ifile = beam['beams'][field]['i_file']
        outdir = f"{os.path.dirname(ifile)}"
        prefix = f"{outdir}/{cname}"

        rm1dfiles = doc['rm1dfiles']
        fdfFile = f"{outdir}/{rm1dfiles['FDF_dirty']}"
        rmsfFile = f"{outdir}/{rm1dfiles['RMSF']}"
        weightFile = f"{outdir}/{rm1dfiles['weights']}"
        rmSynthFile = f"{outdir}/{rm1dfiles['summary_json']}"
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
                                 verbose=clargs.rm_verbose)

        if clargs.database:
            # Load into Mongo
            myquery = {"Component_ID": cname}

            newvalues = {"$set": {f"rmclean1d": True}}
            comp_col.update_one(myquery, newvalues)

            newvalues = {"$set": {f"rmclean_summary": outdict}}
            comp_col.update_one(myquery, newvalues)

    except KeyError:
        print('Failed to load data! RM-CLEAN not applied to component!')
        print(f'Island is {iname}, component is {cname}')
        return


def rmclean3d(args):
    island_id, clargs, outdir, host, field, verbose = args

    with pymongo.MongoClient(host=host) as client:
        mydb = client['spiceracs']  # Create/open database
        isl_col = mydb['islands']  # Create/open collection
        comp_col = mydb['components']  # Create/open collection
        beams_col = mydb['beams']  # Create/open collection

    # Basic querey
    myquery = {"Source_ID": island_id}

    doc = comp_col.find_one(myquery)
    iname = doc['Source_ID']
    prefix = f"{iname}_"
    beam = beams_col.find_one({'Source_ID': iname})
    ifile = beam['beams'][field]['i_file']
    outdir = os.path.dirname(ifile)

    island = isl_col.find_one(myquery)
    rm3dfiles = island["rm3dfiles"]

    cleanFDF, ccArr, iterCountArr, residFDF, headtemp = do_RMclean_3D.run_rmclean(
        fitsFDF=f"{outdir}/{rm3dfiles['FDF_real_dirty']}",
        fitsRMSF=f"{outdir}/{rm3dfiles['RMSF_tot']}",
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
        prefixOut=prefix,
        outDir=outdir,
        write_separate_FDF=True,
        verbose=clargs.rm_verbose)


def main(pool, args, verbose=False):
    outdir = args.outdir
    if outdir[-1] == '/':
        outdir = outdir[:-1]
    outdir = f'{outdir}/cutouts'
    field = args.field

    host = args.host
    # default connection (ie, local)
    with pymongo.MongoClient(host=host) as client:
        mydb = client['spiceracs']  # Create/open database
        isl_col = mydb['islands']  # Create/open collection
        comp_col = mydb['components']  # Create/open collection
        beams_col = mydb['beams']  # Create/open collection

    query = {
        '$and':  [
            {f'beams.{field}': {'$exists': True}},
            {f'beams.{field}.DR1': True}
        ]
    }

    beams = beams_col.find(query).sort('Source_ID')
    all_island_ids = sorted(beams_col.distinct('Source_ID', query))

    query = {
        '$and': [
            {
                'Source_ID': {'$in': all_island_ids}
            },
            {
                'rmsynth3d': True
            }
        ]
    }

    islands = isl_col.find(query).sort('Source_ID')
    island_ids = [doc['Source_ID'] for doc in islands]
    n_island = isl_col.count_documents(query)

    query = {
        '$and': [
            {
                'Source_ID': {'$in': all_island_ids}
            },
            {
                'rmsynth1d': True
            }
        ]
    }

    components = comp_col.find(query).sort('Source_ID')
    component_ids = [doc['Component_ID'] for doc in components]
    n_comp = comp_col.count_documents(query)

    if args.limit is not None:
        count = args.limit
        n_comp = count
        n_island = count
        island_ids = island_ids[:count]
        component_ids = component_ids[:count]

    if args.dimension == '1d':
        if verbose:
            print(f'Running RM-CLEAN on {n_comp} components')
        inputs = [[comp_id, args, outdir, host, field, verbose]
                  for i, comp_id in enumerate(component_ids)]
        if (pool.__class__.__name__ is 'MPIPool' or
                pool.__class__.__name__ is 'SerialPool'):
            if verbose:
                print('Running 1D RM-CLEAN...')
            tic = time.perf_counter()
            list(pool.map(rmclean1d, inputs))
            toc = time.perf_counter()
            if verbose:
                print(f'Time taken was {toc - tic}s')

        elif pool.__class__.__name__ is 'MultiPool':
            list(tqdm(
                pool.imap_unordered(rmclean1d, inputs),
                total=n_comp,
                desc='Running 1D RM-CLEAN',
                disable=(not verbose)
            )
            )

    elif args.dimension == '3d':
        if verbose:
            print(f'Running RM-CLEAN on {n_island} islands')
        inputs = [[island_id, args, outdir, host, field, verbose]
                  for i, island_id in enumerate(island_ids)]
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
                total=n_island,
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
        "field",
        metavar="field",
        type=str,
        help="RACS field to mosaic - e.g. 2132-50A."
    )
    parser.add_argument(
        'outdir',
        metavar='outdir',
        type=str,
        help='Directory containing cutouts (in subdir outdir/cutouts).')

    parser.add_argument(
        'host',
        metavar='host',
        type=str,
        help='Host of mongodb (probably $hostname -i).')

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

    parser.add_argument("--validate", dest="validate", action="store_true",
                        help="Run on Stokes I [False].")

    parser.add_argument("--limit", dest="limit", default=None,
                        type=int, help="Limit number of sources [All].")

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

    host = args.host
    if verbose:
        print('Testing MongoDB connection...')
    # default connection (ie, local)
    with pymongo.MongoClient(host=host) as client:
        try:
            client.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")
        else:
            if verbose:
                print('MongoDB connection succesful!')

    main(pool, args, verbose=verbose)


if __name__ == "__main__":
    cli()
