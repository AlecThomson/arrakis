#!/usr/bin/env python
from utils import getdata, MyEncoder, head2dict
import numpy as np
from tqdm import trange, tqdm
import sys
import os
from dataclasses import dataclass, asdict, make_dataclass
import dataclasses
from astropy.io.fits import Header
import json
import pymongo
from astropy.io import fits
from astropy.table import Table
import time
import warnings
import functools
print = functools.partial(print, flush=True)


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
    if outdir[-1] == '/':
        outdir = outdir[:-1]
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

    # Init cutouts
    i_cutouts = []
    q_cutouts = []
    u_cutouts = []
    if datadict['v_cube'] is not None:
        v_cutouts = []
    else:
        v_cutouts = None

    outdir = f'{outdir}/cutouts'
    if dryrun:
        if verbose:
            print('Dry run -- not saving to disk.')
    else:
        if verbose:
            print(f'Saving to {outdir}/')
    source_dict_list = []

    i_cube = datadict['i_cube']
    q_cube = datadict['q_cube']
    u_cube = datadict['u_cube']
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

    for i in trange(
        count,
        total=count,
        disable=(not verbose),
        desc='Extracting cubelets'
    ):
        if loners:
            if datadict['i_tab']['col_n_components'][i] > 1:
                continue

        # Catch out 0-component sources
        try:
            components = datadict['i_tab_comp'][datadict['i_tab_comp']
                                                ['col_island_id'] == datadict['i_tab']['col_island_id'][i]]

        except KeyError:
            if verbose:
                print(
                    f"No matches found for key {datadict['i_tab']['col_island_id'][i]}")
                print(
                    f"Number of components is {datadict['i_tab']['col_n_components'][i]}!")
            continue
        source_dict = {}
        # Skip if source is outside of cube bounds
        if (y_max[i] > i_cube.shape[2] or x_max[i] > i_cube.shape[2] or
                x_min[i] < 0 or y_min[i] < 0):
            continue

        starty = int(y_min[i]-pad*pixels_per_beam)
        stopy = int(y_max[i]+pad*pixels_per_beam)
        startx = int(x_min[i]-pad*pixels_per_beam)
        stopx = int(x_max[i]+pad*pixels_per_beam)

        # Check if pad puts bbox outside of cube
        if starty < 0:
            starty = y_min[i]
        if startx < 0:
            startx = x_min[i]
        if stopy > i_cube.shape[2]:
            stopy = y_max[i]
        if stopx > i_cube.shape[2]:
            stopx = x_max[i]

        i_cutout = i_cube[:, starty:stopy, startx:stopx]
        q_cutout = q_cube[:, starty:stopy, startx:stopx]
        u_cutout = u_cube[:, starty:stopy, startx:stopx]

        i_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
        q_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
        u_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]

        i_cutouts.append(i_cutout)
        q_cutouts.append(q_cutout)
        u_cutouts.append(u_cutout)

        if v_cube is not None:
            v_cutout = v_cube[:, starty:stopy, startx:stopx]
            v_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
            v_cutouts.append(v_cutout)

        source_dict['header'] = i_cutout.header
        del source_dict['header']['SLICE']
        for name in datadict['i_tab'].columns:
            source_dict[name.replace('col_', '')
                        ] = datadict['i_tab'][name][i]

        for comp, row in enumerate(components.iterrows()):
            source_dict[f'component_{comp+1}'] = {}
            for name in row[1].keys():
                if type(row[1][name]) is bytes:
                    source_dict[f'component_{comp+1}'][name.replace('col_', '')
                                                       ] = row[1][name]
                else:
                    source_dict[f'component_{comp+1}'][name.replace('col_', '')
                                                       ] = row[1][name]

        source_dict_list.append(source_dict)

    # Set up locations where files will be saved
    for i in trange(
        len(source_dict_list),
        disable=(not verbose),
        desc='Finding locations'
    ):
        for stoke in ['i', 'q', 'u', 'v']:
            if stoke == 'v' and datadict['v_cube'] is None:
                pass
            else:
                name = source_dict_list[i]['island_name']
                outname = f'{name}.cutout.{stoke}.fits'
                source_dict_list[i][f'{stoke}_file'] = outname

    if v_cube is not None:
        stoke = 'v'
        name = source_dict_list[i]['island_name']
        outname = f'{name}.cutout.{stoke}.fits'
        source_dict_list[i][f'{stoke}_file'] = outname

    cutouts = {
        "i": i_cutouts,
        "q": q_cutouts,
        "u": u_cutouts
    }

    if v_cube is not None:
        cutouts.update(
            {
                "v": v_cutouts
            }
        )

    return cutouts, source_dict_list, outdir


def getbytes(cutout):
    """Worker task for size estimate.

    Args:
        cutout (SpectralCube): The cutout around a source.

    Returns:
        mbytes (float): The size of the cutout in MB

    """

    mbytes = cutout[0, :, :].nbytes*cutout.shape[0]*1e-6
    return mbytes


def getsize(pool, cutouts, verbose=True):
    """Find the size of all cutouts in a cube.

    Args:
        pool (Pool): The Shwimmbad or Multiprocessing pool for
            parallelisation.
        cutouts (dict): Contains lists of I, Q, and U cutouts.

    Kwargs:
        verbose (bool): Print out messages.

    """
    if (pool.__class__.__name__ is 'MPIPool' or
            pool.__class__.__name__ is 'SerialPool'):
        if verbose:
            print(f'Getting sizes...')
        tic = time.perf_counter()
        sizes_bytes = list(
            pool.map(getbytes, [cutouts[i] for i in range(len(cutouts))]
                     )
        )
        toc = time.perf_counter()
        if verbose:
            print(f'Time taken was {toc - tic}s')

    elif pool.__class__.__name__ is 'MultiPool':
        sizes_bytes = list(tqdm(
            pool.imap_unordered(getbytes, [cutouts[i]
                                           for i in range(len(cutouts))]),
            total=len(cutouts),
            desc='Getting sizes',
            disable=(not verbose)
        )
        )

    sizes_bytes = np.array(sizes_bytes)
    print('Size in MB: ', sizes_bytes.sum())
    print('Size in GB: ', sizes_bytes.sum()/1000)


def writefits(arg):
    """Worker that cutouts to disk.

    Writes a cutout, as stored in source_dict, to disk. The file
    location should already be specified in source_dict. This format is
    intended for parallel use with pool.map syntax.

    Args:
        arg: The tuple of (source_dict, cutout, stoke)
            source_dict (dict): Source metadata.
            cutout (SpectralCube): Cutout data to write.
            stoke (str): Which Stokes to write. Is a string of either 'i',
                'q', 'u', or 'v'.
            outdir (str): Location of output directory.

    """

    outfile, idx, plane = arg

    # Reopen file and update
    outfh = fits.open(outfile, mode='update')
    outfh[0].data[idx] = plane
    outfh.flush()


def writeloop(pool, cutouts, source_dict_list, datadict, outdir, verbose=True):
    """Main loop for writing to disk

    Loops over all sources in I, Q, and U, and writes to disk in
    parallel.

    Args:
        pool (Pool): The Shwimmbad or Multiprocessing pool for
            parallelisation.
        cutouts (dict): Contains lists of I, Q, and U cutouts.
        source_dict_list (list): List of dictionaries which contain the
            metadata of each source.
        datadict (dict): Dictionary containing tables and cubes.
        outdir (str): Location of output directory.

    Kwargs:
        verbose (bool): Print out messages.

    """
    for stoke in ['i', 'q', 'u', 'v']:
        if stoke == 'v' and datadict['v_cube'] is None:
            pass
        else:
            for i, cutout in enumerate(tqdm(cutouts[stoke], desc=f'Stokes {stoke}', disable=(not verbose))):
                outfile = f"{outdir}/{source_dict_list[i][f'{stoke}_file']}"
                # Write blank file
                blank = np.zeros(cutout.shape)
                hdu = fits.PrimaryHDU(blank)
                hdulist = fits.HDUList([hdu])
                hdulist.writeto(outfile, overwrite=True)
                inputs = [[outfile, idx, cutout.unmasked_data[idx]]
                          for idx in range(len(cutout))]

                if (pool.__class__.__name__ is 'MPIPool' or
                        pool.__class__.__name__ is 'SerialPool'):
                    if verbose:
                        print(f'Writing Stokes {stoke}...')
                    tic = time.perf_counter()
                    list(pool.map(writefits, inputs))
                    toc = time.perf_counter()
                    if verbose:
                        print(f'Time taken was {toc - tic}s')

                elif pool.__class__.__name__ is 'MultiPool':
                    list(tqdm(pool.imap_unordered(writefits, inputs), total=len(
                        cutout.spectral_axis),desc=f'Writing per channel',  disable=(not verbose)))

                outfh = fits.open(outfile, mode='update')
                outfh[0].header = cutout.header
                del outfh[0].header['SLICE']
                outfh.flush()




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
    cutouts, source_dict_list, outdir = makecutout(pool,
                                                   datadict,
                                                   outdir=outdir,
                                                   pad=pad,
                                                   dryrun=dryrun,
                                                   limit=args.limit,
                                                   loners=args.loners,
                                                   verbose=verbose
                                                   )

    # Check size of cube
    if args.getsize:
        if verbose:
            print('Checking size of single cube...')
        getsize(pool, cutouts['i'], verbose=verbose)

    # Write to disk
    if not dryrun:
        if verbose:
            print('Writing to disk...')
        writeloop(pool, cutouts, source_dict_list,
                  datadict, outdir, verbose=verbose)

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
