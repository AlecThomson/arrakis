#!/usr/bin/env python
from spiceracs.utils import getdata
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
import time
import warnings
import functools
print = functools.partial(print, flush=True)



def makecutout(pool, datadict, outdir='.', pad=0, dryrun=False, verbose=True):
    """Main cutout task.

    Takes in table data from Selavy, and cuts out data cubes using
    SpectralCube.

    Args:
        datadict (dict): Dictionary containing tables and cubes.

    Kwargs:
        outdir (str): Directory to save data to.
        pad (float): Fractional padding around Selavy islands to cut out.
        dryrun (bool): Whether to do a dry-run (save no files).
        verbose (bool): Print out messages.

    Returns:
        cutouts (dict): Contains lists of I, Q, and U cutouts.
        source_dict_list (list): List of dictionaries which contain the
            metadata of each source.

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

    # Init cutouts
    i_cutouts = []
    q_cutouts = []
    u_cutouts = []
    outdir = f'{outdir}/cutouts'
    if dryrun:
        if verbose: print('Dry run -- not saving to disk.')
    else:
        if verbose: print(f'Saving to {outdir}/')
    source_dict_list = []

    i_cube = datadict['i_cube']
    q_cube = datadict['q_cube']
    u_cube = datadict['u_cube']
    if not dryrun:
        try:
            os.mkdir(outdir)
            print('Made directory.')
        except FileExistsError:
            print('Directory exists.')

    # TO-DO: Max cut on size
    for i in trange(
        len(datadict['i_tab']),
        total = len(datadict['i_tab']),
        disable=(not verbose),
        desc='Extracting cubelets'
        ):
        source_dict = {}
        # Skip if source is outside of cube bounds
        if (y_max[i] > i_cube.shape[2] or x_max[i] > i_cube.shape[2] or
                x_min[i] < 0 or y_min[i] < 0):
            continue

        # Check if pad puts bbox outside of cube
        elif (int(y_min[i]-pad*dy[i]) > 0 and
              int(x_min[i]-pad*dx[i]) > 0 and
              int(y_max[i]+pad*dy[i]) < i_cube.shape[2] and
              int(x_max[i]+pad*dx[i]) < i_cube.shape[2]):

            i_cutout = i_cube[
                :,
                int(y_min[i]-pad*dy[i]):int(y_max[i]+pad*dy[i]),
                int(x_min[i]-pad*dx[i]):int(x_max[i]+pad*dx[i])
            ]
            q_cutout = q_cube[
                :,
                int(y_min[i]-pad*dy[i]):int(y_max[i]+pad*dy[i]),
                int(x_min[i]-pad*dx[i]):int(x_max[i]+pad*dx[i])
            ]
            u_cutout = u_cube[
                :,
                int(y_min[i]-pad*dy[i]):int(y_max[i]+pad*dy[i]),
                int(x_min[i]-pad*dx[i]):int(x_max[i]+pad*dx[i])
            ]

            i_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
            q_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
            u_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
            i_cutouts.append(i_cutout)
            q_cutouts.append(q_cutout)
            u_cutouts.append(u_cutout)

            source_dict['header'] = i_cutout.header
            for name in datadict['i_tab'].colnames:
                source_dict[name.replace('col_', '')
                            ] = datadict['i_tab'][name][i]
            source_dict_list.append(source_dict)
        else:
            i_cutout = i_cube[:, y_min[i]:y_max[i], x_min[i]:x_max[i]]
            q_cutout = q_cube[:, y_min[i]:y_max[i], x_min[i]:x_max[i]]
            u_cutout = u_cube[:, y_min[i]:y_max[i], x_min[i]:x_max[i]]
            i_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
            q_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
            u_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
            i_cutouts.append(i_cutout)
            q_cutouts.append(q_cutout)
            u_cutouts.append(u_cutout)

            source_dict['header'] = i_cutout.header
            for name in datadict['i_tab'].colnames:
                source_dict[name.replace('col_', '')
                            ] = datadict['i_tab'][name][i]
            source_dict_list.append(source_dict)

    # Set up locations where files will be saved
    for i in trange(
        len(source_dict_list),
        disable=(not verbose),
        desc='Finding locations'
    ):
        for stoke in ['i', 'q', 'u']:
            name = source_dict_list[i]['island_name']
            outname = f'{outdir}/{name}.cutout.{stoke}.fits'
            source_dict_list[i][f'{stoke}_file'] = outname

    cutouts = {
        "i": i_cutouts,
        "q": q_cutouts,
        "u": u_cutouts
    }

    return cutouts, source_dict_list


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
        if verbose: print(f'Getting sizes...')
        tic = time.perf_counter()
        sizes_bytes = list(
            pool.map(getbytes, [cutouts[i] for i in range(len(cutouts))]
                     )
        )
        toc = time.perf_counter()
        if verbose: print(f'Time taken was {toc - tic}s')

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
                'q', or 'u'.

    """
    source_dict, cutout, stoke = arg
    outfile = source_dict[f'{stoke}_file']
    cutout.write(outfile, format='fits', overwrite=True)


def writeloop(pool, cutouts, source_dict_list, verbose=True):
    """Main loop for writing to disk

    Loops over all sources in I, Q, and U, and writes to disk in
    parallel.

    Args:
        pool (Pool): The Shwimmbad or Multiprocessing pool for
            parallelisation.
        cutouts (dict): Contains lists of I, Q, and U cutouts.
        source_dict_list (list): List of dictionaries which contain the
            metadata of each source.

    Kwargs:
        verbose (bool): Print out messages.

    """
    for stoke in ['i', 'q', 'u']:
        if (pool.__class__.__name__ is 'MPIPool' or
            pool.__class__.__name__ is 'SerialPool'):
            if verbose: print(f'Writing Stokes {stoke}...')
            tic = time.perf_counter()
            pool.map(
                writefits,
                [[source_dict_list[i], cutouts[stoke][i], stoke]
                    for i in range(len(cutouts[stoke]))]
            )
            toc = time.perf_counter()
            if verbose: print(f'Time taken was {toc - tic}s')

        elif pool.__class__.__name__ is 'MultiPool':
            list(tqdm(
                pool.imap_unordered(
                    writefits,
                    [[source_dict_list[i], cutouts[stoke][i], stoke]
                     for i in range(len(cutouts[stoke]))]
                ),
                total=len(cutouts[stoke]),
                desc=f'Stokes {stoke}',
                disable=(not verbose)
            )
            )


def head2dict(h):
    """Convert FITS header to a dict.

    Writes a cutout, as stored in source_dict, to disk. The file location
    should already be specified in source_dict. This format is intended
    for parallel use with pool.map syntax.

    Args:
        h: An astropy FITS header.

    Returns:
        data (dict): The FITS head converted to a dict.

    """
    data = {}
    for c in h.__dict__['_cards']:
        if c[0] == '':
            continue
        data[c[0]] = c[1]
    return(data)


class MyEncoder(json.JSONEncoder):
    """Cutom JSON encorder.

    Parses the data stored in source_dict to JSON without
    errors.

    """

    def default(self, obj): # pylint: disable=E0202
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, fits.Header):
            return head2dict(obj)
        elif dataclasses.is_dataclass(obj):
            return dataclasses.asdict(obj)
        else:
            return super(MyEncoder, self).default(obj)


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
    if verbose: print('Total documents:', count)
    client.close()


def main(pool, args, verbose=True):
    """Main script.

    """
    # Sort out args
    cubedir = args.cubedir
    tabledir = args.tabledir
    if cubedir[-1] == '/':
        cubedir = cubedir[:-1]

    if tabledir[-1] == '/':
        tabledir = tabledir[:-1]

    # Read in data
    if verbose: print('Reading data...')
    datadict = getdata(cubedir, tabledir, verbose=verbose)

    # Make cutouts
    outdir = args.outdir
    if outdir[-1] == '/':
        outdir = outdir[:-1]
    pad = args.pad
    dryrun = args.dryrun
    if verbose: print('Making cutouts....')
    cutouts, source_dict_list = makecutout(pool,
        datadict,
        outdir=outdir,
        pad=pad,
        dryrun=dryrun,
        verbose=verbose
    )

    # Check size of cube
    if args.getsize:
        if verbose: print('Checking size of single cube...')
        getsize(pool, cutouts['i'], verbose=verbose)

    # Write to disk
    if not dryrun:
        if verbose: print('Writing to disk...')
        writeloop(pool, cutouts, source_dict_list, verbose=verbose)

    # Update MongoDB
    if args.database:
        if verbose: print('Updating MongoDB...')
        database(source_dict_list, verbose=True)

    if verbose: print('Done!')

def cli():
    """Command-line interface.
    """
    import argparse
    import schwimmbad
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
    SPICE-RACS Stage 1:
    Produce cubelets from a RACS field using a Selavy table.

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
        'outdir',
        metavar='outdir',
        type=str,
        help='Directory to store cutouts.')

    parser.add_argument(
        'pad',
        metavar='pad',
        type=float,
        default=0,
        help='Fractional padding around islands [0 -- no padding].')

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
    pool =  schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)
    if args.mpi:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
    if verbose:
        print(f"Using pool: {pool.__class__.__name__}")

    if args.database:
        if verbose: print('Testing MongoDB connection...')
        client = pymongo.MongoClient()  # default connection (ie, local)
        try:
            client.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")
        else:
            if verbose: print('MongoDB connection succesful!')
        client.close()

    main(pool, args, verbose=verbose)
    pool.close()

if __name__ == "__main__":
    cli()