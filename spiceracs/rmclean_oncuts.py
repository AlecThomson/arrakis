#!/usr/bin/env python3
from spiceracs.utils import getfreq, MyEncoder, tqdm_dask
import json
import numpy as np
import os
import pymongo
import sys
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
from IPython import embed
import dask
from dask import delayed
from dask.distributed import Client, progress, LocalCluster
from dask.diagnostics import ProgressBar


@delayed
def rmclean1d(comp_id,
              outdir,
              host,
              field,
              cutoff=-3,
              maxIter=10000,
              gain=0.1,
              showPlots=False,
              database=False,
              rm_verbose=True,
              ):
    """1D RM-CLEAN

    Args:
        comp_id (str): RACS component ID
        host (str): MongoDB host
        field (str): RACS field
        cutoff (int, optional): CLEAN cutoff. Defaults to -3.
        maxIter (int, optional): CLEAN max iterations. Defaults to 10000.
        gain (float, optional): CLEAN gain. Defaults to 0.1.
        showPlots (bool, optional): Show plots. Defaults to False.
        database (bool, optional): Update MongoDB. Defaults to False.
        rm_verbose (bool, optional): Verbose RM-CLEAN. Defaults to True.
    """
    with pymongo.MongoClient(host=host, connect=False) as dbclient:
        mydb = dbclient['spiceracs']  # Create/open database
        isl_col = mydb['islands']  # Create/open collection
        comp_col = mydb['components']  # Create/open collection
        beams_col = mydb['beams']  # Create/open collection

    # Basic querey
    myquery = {"Gaussian_ID": comp_id}

    doc = comp_col.find_one(myquery)

    iname = doc['Source_ID']
    beam = beams_col.find_one({'Source_ID': iname})

    if rm_verbose:
        print(f'Working on {comp_id}')
    try:
        cname = comp_id

        rm1dfiles = doc['rm1dfiles']
        fdfFile = os.path.join(outdir, f"{rm1dfiles['FDF_dirty']}")
        rmsfFile = os.path.join(outdir, f"{rm1dfiles['RMSF']}")
        weightFile = os.path.join(outdir, f"{rm1dfiles['weights']}")
        rmSynthFile = os.path.join(outdir, f"{rm1dfiles['summary_json']}")

        prefix = os.path.join(os.path.abspath(os.path.dirname(fdfFile)), cname)

        # Sanity checks
        for f in [weightFile, fdfFile, rmsfFile, rmSynthFile]:
            if not os.path.exists(f):
                print("File does not exist: '{:}'.".format(f), end=' ')
                sys.exit()
        nBits = 32
        mDict, aDict = do_RMclean_1D.readFiles(
            fdfFile, rmsfFile, weightFile, rmSynthFile, nBits)

        # Run RM-CLEAN on the spectrum
        outdict, arrdict = do_RMclean_1D.run_rmclean(mDict=mDict,
                                                     aDict=aDict,
                                                     cutoff=cutoff,
                                                     maxIter=maxIter,
                                                     gain=gain,
                                                     nBits=nBits,
                                                     showPlots=showPlots,
                                                     verbose=rm_verbose)

        # Save output
        do_RMclean_1D.saveOutput(outdict,
                                 arrdict,
                                 prefixOut=prefix,
                                 verbose=rm_verbose)

        if database:
            # Load into Mongo
            myquery = {"Gaussian_ID": cname}

            newvalues = {"$set": {f"rmclean1d": True}}
            comp_col.update_one(myquery, newvalues)

            newvalues = {"$set": {f"rmclean_summary": outdict}}
            comp_col.update_one(myquery, newvalues)

    except KeyError:
        print('Failed to load data! RM-CLEAN not applied to component!')
        print(f'Island is {iname}, component is {cname}')
        return


@delayed
def rmclean3d(island_id,
              outdir,
              host,
              cutoff=-3,
              maxIter=10000,
              gain=0.1,
              database=False,
              rm_verbose=False,
              ):
    """3D RM-CLEAN

    Args:
        island_id (str): RACS Island ID
        host (str): MongoDB host
        field (str): RACS field
        cutoff (int, optional): CLEAN cutoff. Defaults to -3.
        maxIter (int, optional): CLEAN max iterations. Defaults to 10000.
        gain (float, optional): CLEAN gain. Defaults to 0.1.
        rm_verbose (bool, optional): Verbose RM-CLEAN. Defaults to False.
    """
    with pymongo.MongoClient(host=host, connect=False) as dbclient:
        mydb = dbclient['spiceracs']  # Create/open database
        isl_col = mydb['islands']  # Create/open collection
        comp_col = mydb['components']  # Create/open collection

    # Basic querey
    myquery = {"Source_ID": island_id}

    doc = comp_col.find_one(myquery)
    iname = doc['Source_ID']
    prefix = f"{iname}_"

    island = isl_col.find_one(myquery)
    rm3dfiles = island["rm3dfiles"]

    cleanFDF, ccArr, iterCountArr, residFDF, headtemp = do_RMclean_3D.run_rmclean(
        fitsFDF=os.path.join(outdir, rm3dfiles['FDF_real_dirty']),
        fitsRMSF=os.path.join(outdir, rm3dfiles['RMSF_tot']),
        cutoff=cutoff,
        maxIter=maxIter,
        gain=gain,
        chunksize=None,
        nBits=32,
        verbose=rm_verbose)

    # Write results to disk
    do_RMclean_3D.writefits(
        cleanFDF,
        ccArr,
        iterCountArr,
        residFDF,
        headtemp,
        prefixOut=prefix,
        outDir=os.path.abspath(os.path.dirname(
            os.path.join(outdir, rm3dfiles['FDF_real_dirty']))),
        write_separate_FDF=True,
        verbose=rm_verbose)
    
    if database:
        # Load into Mongo
        myquery = {"Source_ID": iname}
        newvalues = {"$set": {f"rmclean3d": True}}
        isl_col.update_one(myquery, newvalues)


def main(field,
         outdir,
         host,
         client,
         dimension='1d',
         verbose=True,
         database=False,
         validate=False,
         limit=None,
         cutoff=-3,
         maxIter=10000,
         gain=0.1,
         showPlots=False,
         rm_verbose=False
         ):
    """Main script

    Args:
        field (str): RACS field
        outdir (str): Work directory (contains 'cutouts' as subdir)
        host (str): MongoDB host
        client (Client): Dask client
        dimension (str, optional): RM-CLEAN dimension. Defaults to '1d'.
        verbose (bool, optional): Verbose output. Defaults to True.
        database (bool, optional): Update MongoDB. Defaults to False.
        validate (bool, optional): Run on Stokes I. Defaults to False.
        limit (int, optional): Limit number of sources to CLEAN. Defaults to None.
        cutoff (float, optional): CLEAN cutof. Defaults to -3.
        maxIter (int, optional): CLEAN max iterations. Defaults to 10000.
        gain (float, optional): CLEAN gain. Defaults to 0.1.
        showPlots (bool, optional): Show CLEAN plots. Defaults to False.
        rm_verbose (bool, optional): Verbose RM-CLEAN. Defaults to False.
    """
    if outdir[-1] == '/':
        outdir = outdir[:-1]
    outdir = f'{outdir}/cutouts'

    # default connection (ie, local)
    with pymongo.MongoClient(host=host, connect=False) as dbclient:
        mydb = dbclient['spiceracs']  # Create/open database
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
    component_ids = [doc['Gaussian_ID'] for doc in components]
    n_comp = comp_col.count_documents(query)

    if limit is not None:
        count = limit
        n_comp = count
        n_island = count
        island_ids = island_ids[:count]
        component_ids = component_ids[:count]

    outputs = []
    if dimension == '1d':
        if verbose:
            print(f'Running RM-CLEAN on {n_comp} components')
        for i, comp_id in enumerate(component_ids):
            if i > n_comp+1:
                break
            else:
                output = rmclean1d(comp_id,
                                   outdir,
                                   host,
                                   field,
                                   cutoff=cutoff,
                                   maxIter=maxIter,
                                   gain=gain,
                                   showPlots=showPlots,
                                   database=database,
                                   rm_verbose=rm_verbose)
                outputs.append(output)

    elif dimension == '3d':
        if verbose:
            print(f'Running RM-CLEAN on {n_island} islands')

        for i, island_id in enumerate(island_ids):
            if i > n_island+1:
                break
            else:
                output = rmclean3d(island_id=island_id,
                                   outdir=outdir,
                                   host=host,
                                   cutoff=cutoff,
                                   maxIter=maxIter,
                                   gain=gain,
                                   database=database,
                                   rm_verbose=rm_verbose)
                outputs.append(output)

    results = client.persist(outputs)
    tqdm_dask(results, desc='Running RM-CLEAN', disable=(not verbose))

    if verbose:
        print('Done!')


def cli():
    """Command-line interface
    """
    import argparse
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
    parser.add_argument("-n", dest="maxIter", type=int, default=10000,
                        help="maximum number of CLEAN iterations [10000].")
    parser.add_argument("-g", dest="gain", type=float, default=0.1,
                        help="CLEAN loop gain [0.1].")
    parser.add_argument("-p", dest="showPlots", action="store_true",
                        help="show the plots [False].")
    parser.add_argument("-rmv", dest="rm_verbose", action="store_true",
                        help="Verbose RM-CLEAN [False].")

    args = parser.parse_args()

    cluster = LocalCluster(n_workers=20)
    client = Client(cluster)

    verbose = args.verbose
    host = args.host
    if verbose:
        print('Testing MongoDB connection...')
    # default connection (ie, local)
    with pymongo.MongoClient(host=host, connect=False) as dbclient:
        try:
            dbclient.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")
        else:
            if verbose:
                print('MongoDB connection succesful!')

    main(field=args.field,
         outdir=args.outdir,
         host=host,
         client=client,
         dimension=args.dimension,
         verbose=verbose,
         database=args.database,
         validate=args.validate,
         limit=args.limit,
         cutoff=args.cutoff,
         maxIter=args.maxIter,
         gain=args.gain,
         showPlots=args.showPlots,
         rm_verbose=args.rm_verbose,
         )


if __name__ == "__main__":
    cli()
