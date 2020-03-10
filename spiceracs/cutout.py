#!/usr/bin/env python
from utils import getdata, MyEncoder
import numpy as np
from tqdm import trange, tqdm
import sys
import os
import subprocess
import shlex
import json
import pymongo
from astropy.io import fits
import time
import warnings
import functools
print = functools.partial(print, flush=True)


def cutout_worker(args):
    i, i_tab, i_tab_comp, i_file, x_min, x_max, y_min, y_max, pad, pixels_per_beam, shape, outdir, v_cube_is_None, dryrun, loners, verbose = args
    if loners:
        if i_tab['col_n_components'][i] > 1:
            return

    # Catch out 0-component sources
    try:
        components = i_tab_comp[i_tab_comp['col_island_id']
                                == i_tab['col_island_id'][i]]

    except KeyError:
        if verbose:
            print(
                f"No matches found for key {i_tab['col_island_id'][i]}")
            print(
                f"Number of components is {i_tab['col_n_components'][i]}!")
        return
    source_dict = {}
    # Skip if source is outside of cube bounds
    if (y_max[i] > shape[2] or x_max[i] > shape[2] or
            x_min[i] < 0 or y_min[i] < 0):
        return

    starty = int(y_min[i]-pad*pixels_per_beam)
    stopy = int(y_max[i]+pad*pixels_per_beam)
    startx = int(x_min[i]-pad*pixels_per_beam)
    stopx = int(x_max[i]+pad*pixels_per_beam)

    # Check if pad puts bbox outside of cube
    if starty < 0:
        starty = y_min[i]
    if startx < 0:
        startx = x_min[i]
    if stopy > shape[2]:
        stopy = y_max[i]
    if stopx > shape[2]:
        stopx = x_max[i]

    for comp, row in enumerate(components.iterrows()):
        source_dict[f'component_{comp+1}'] = {}
        for name in row[1].keys():
            if type(row[1][name]) is bytes:
                source_dict[f'component_{comp+1}'][name.replace('col_', '')
                                                   ] = row[1][name]
            else:
                source_dict[f'component_{comp+1}'][name.replace('col_', '')
                                                   ] = row[1][name]

    for name in i_tab.columns:
        source_dict[name.replace('col_', '')
                    ] = i_tab[name][i]

    for stoke in ['i', 'q', 'u', 'v']:
        if stoke == 'v' and v_cube_is_None:
            pass
        else:
            name = source_dict['island_name']
            # if stoke == 'i':
            #    print(f"I'm working on Island {name}")
            outname = f'test.{name}.cutout.{stoke}.fits'
            source_dict[f'{stoke}_file'] = outname
        outfile = f"{outdir}/{outname}"
        command = f"fitscopy {i_file}[{startx+1}:{stopx},{starty+1}:{stopy}] {outfile}"
        command = shlex.split(command)
        if not dryrun:
            proc = subprocess.run(command, capture_output=(
                not verbose), encoding="utf-8", check=True)
            headfile = f'{outdir}/{name}.cutout.i.fits'
            source_dict['header'] = fits.getheader(headfile)
            del source_dict['header']['SLICE']

    return source_dict


def makecutout(pool, datadict, outdir='.', pad=0, dryrun=False, limit=None, loners=False, verbose=True):
    """Main cutout task.

    Takes in table data from Selavy, and cuts out data cubes using
    SpectralCube.

    Args:
        datadict (dict): Dictionary containing tables and cubes.

    Kwargs:
        outdir (str): Directory to save data to.
        pad (float): Fractional padding around Selavy islands to cut out.
        dryrun (bool): Whether to do a dry-run (save no files).
        loners (bool): Only run on single-component islands.
        limit (int): Number of sources to cut out.
        verbose (bool): Print out messages.

    Returns:
        cutouts (dict): Contains lists of I, Q, and U cutouts.
        source_dict_list (list): List of dictionaries which contain the
            metadata of each source.
        outdir (str): Location of output directory.

    """
    # Get bounding boxes in WCS
    ra_min, dec_min, freq = datadict['wcs_taylor'].all_pix2world(
        datadict['i_tab']['col_x_min'], datadict['i_tab']['col_y_min'], 0, 0)
    ra_max, dec_max, freq = datadict['wcs_taylor'].all_pix2world(
        datadict['i_tab']['col_x_max'], datadict['i_tab']['col_y_max'], 0, 0)

    # Get bounding boxes in cube pixels
    x_min, y_min, _ = np.array(datadict['wcs_cube'].all_world2pix(
        ra_min, dec_min, freq, 0)).astype(int)
    x_max, y_max, _ = np.array(datadict['wcs_cube'].all_world2pix(
        ra_max, dec_max, freq, 0)).astype(int)
    dy, dx = y_max - y_min, x_max-x_min

    # Get beam info - use major axis of beam
    pixels_per_beam = int(
        datadict['i_cube'].header['BMAJ']/datadict['i_cube'].header['CDELT2'])

    outdir = f'{outdir}/cutouts'
    if dryrun:
        if verbose:
            print('Dry run -- not saving to disk.')
    else:
        if verbose:
            print(f'Saving to {outdir}/')
    if datadict['v_cube'] is not None:
        v_cube = datadict['v_cube']
    else:
        v_cube = None

    if not dryrun:
        try:
            os.mkdir(outdir)
            print('Made directory.')
        except FileExistsError:
            print('Directory exists.')

    if limit is not None:
        count = limit
        if verbose:
            print(f'Only cutting out {limit} sources...')
    else:
        count = len(datadict['i_tab'])

    if loners and verbose:
        sub_count = len(datadict['i_tab'][:count]['col_n_components']
                        [datadict['i_tab'][:count]['col_n_components'] < 2])
        print(
            f'Only cutting out {sub_count} sources (skipping muti-component sources)..')

    # datadict['i_tab_comp'].add_index('col_island_id')

    # for i in trange(
    #    count,
    #    total=count,
    #    disable=(not verbose),
    #    desc='Extracting cubelets'
    # ):

    inputs = [[i, datadict['i_tab'], datadict['i_tab_comp'], datadict['i_file'],
               x_min, x_max, y_min, y_max, pad, pixels_per_beam,
               datadict['i_cube'].shape, outdir, (v_cube is None), dryrun, loners, verbose] for i in range(count)]

    if pool.__class__.__name__ is 'MPIPool':
        if verbose:
            print('Extracting cubelets...')
        tic = time.perf_counter()
        source_dict_list = list(pool.map(cutout_worker, inputs))
        toc = time.perf_counter()
        if verbose:
            print(f'Time taken was {toc - tic}s')

    elif pool.__class__.__name__ is 'SerialPool':
        source_dict_list = []
        for i in trange(count):
            source_dict_list.append(cutout_worker([i, datadict['i_tab'], datadict['i_tab_comp'], datadict['i_file'],
                                                   x_min, x_max, y_min, y_max, pad, pixels_per_beam,
                                                   datadict['i_cube'].shape, outdir, (v_cube is None), dryrun, loners, verbose]))

    elif pool.__class__.__name__ is 'MultiPool':
        source_dict_list = list(tqdm(
            pool.imap_unordered(cutout_worker, inputs),
            total=count,
            desc='Extracting cubelets',
            disable=(not verbose)
        )
        )

    #source_dict_list = list(pool.map(cutout_worker, inputs))

    return source_dict_list


def database(source_dict_list, verbose=True):
    """Add data to MongoDB.

    Args:
        source_dict_list (list): List of dictionaries which contain the
            metadata of each source.

    Kwargs:
        verbose (bool): Print out messages.

    """
    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

    for source_dict in tqdm(
        source_dict_list,
        desc='Loading into DB',
        disable=(not verbose)
    ):
        json_data = json.loads(json.dumps(source_dict, cls=MyEncoder))
        name = source_dict['island_name']
        count = mycol.count_documents({'island_name': name})

        if count > 0:
            # Replace existing
            cursor = mycol.find({'island_name': name})
            for doc in cursor:
                mycol.replace_one({'_id': doc['_id']}, json_data)
        else:
            # Insert new
            mycol.insert_one(json_data)
    # Check if all cutouts are in collection
    count = mycol.count_documents({})
    if verbose:
        print('Total documents:', count)
    client.close()


def main(pool, args, verbose=True):
    """Main script.

    """
    # Sort out args
    cubedir = args.cubedir
    tabledir = args.tabledir
    mapdata = args.mapdata

    # Read in data
    if verbose:
        print('Reading data...')
    datadict = getdata(cubedir, tabledir, mapdata, verbose=verbose)

    # Make cutouts
    pad = args.pad
    dryrun = args.dryrun
    if verbose:
        print('Making cutouts....')
    outdir = args.outdir
    if outdir[-1] == '/':
        outdir = outdir[:-1]
    source_dict_list = makecutout(pool,
                                  datadict,
                                  outdir=outdir,
                                  pad=pad,
                                  dryrun=dryrun,
                                  limit=args.limit,
                                  loners=args.loners,
                                  verbose=verbose
                                  )

    # Update MongoDB
    if args.database:
        if verbose:
            print('Updating MongoDB...')
        database(source_dict_list, verbose=True)

    pool.close()

    if verbose:
        print('Done!')


def cli():
    """Command-line interface.
    """
    import argparse
    import schwimmbad
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
    warnings.filterwarnings(
        "ignore", message="Degrees of freedom <= 0 for slice.")

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
    SPICE-RACS Stage 1:
    Produce cubelets from a RACS field using a Selavy table.
    If Stokes V is present, it will be squished into RMS spectra.

    To use with MPI:
       mpirun -n $NPROCS python -u cutout.py $cubedir $tabledir
       $outdir --mpi
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
        'mapdata',
        metavar='mapdata',
        type=str,
        help='2D FITS image corresponding to Selavy table.')

    parser.add_argument(
        'outdir',
        metavar='outdir',
        type=str,
        help='Directory to store cutouts (in subdir outdir/cutouts).')

    parser.add_argument(
        'pad',
        metavar='pad',
        type=float,
        default=0,
        help='Number of beamwidths to pad around source [0 -- no padding].')

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
        "-m",
        dest="database",
        action="store_true",
        help="Add data to MongoDB [False]."
    )

    parser.add_argument("--limit", dest="limit", default=None,
                        type=int, help="Limit number of sources [All].")

    parser.add_argument("--loners", dest="loners", action="store_true",
                        help="Run on single component sources [False].")

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        "--ncores",
        dest="n_cores",
        default=1,
        type=int, help="Number of processes (uses multiprocessing)."
    )
    group.add_argument(
        "--mpi",
        dest="mpi",
        default=False,
        action="store_true",
        help="Run with MPI."
    )

    args = parser.parse_args()
    verbose = args.verbose
    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)
    if args.mpi:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
    print('MPI arg is', args.mpi)
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
