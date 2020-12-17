#!/usr/bin/env python
import warnings
from mpi4py import MPI
from IPython import embed
from functools import partial
import functools
import psutil
import shlex
import subprocess
import schwimmbad
import pymongo
import json
import time
from tqdm import tqdm, trange
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
import sys
import os
from glob import glob
from astropy.table import Table, vstack
from spiceracs.utils import getdata, MyEncoder
from astropy.utils import iers
from spectral_cube import SpectralCube
iers.conf.auto_download = False
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print = functools.partial(
    print, f'[{psutil.Process().cpu_num()},{rank}]', flush=True)
warnings.filterwarnings(
    "ignore", message="Cube is a Stokes cube, returning spectral cube for I component")


def cutout(image, src_name, ra, dec, src_width, outdir, pad=3, verbose=False):
    """Cutout a source from a given image.

    Arguments:
        image {str} -- Name of the FITS image to cutout from
        src_name {str} -- Name of the source
        ra {float} -- RA of source in DEG
        dec {float} -- DEC of source in DEG
        src_width {float} -- Width of source in DEG
        outdir {str} -- Directory to save cutout

    Keyword Arguments:
        pad {int} -- Number of beamwidth to pad cutout (default: {3})
        verbose {bool} -- Verbose output (default: {False})
        dryrun {bool} -- Do everything except the cutout (default: {True})
    """
    if verbose:
        print(f'Reading {image}')
    if outdir[-1] == '/':
        outdir = outdir[:-1]
    basename = os.path.basename(image)
    outname = f'{src_name}.cutout.{basename}'
    outfile = f"{outdir}/{outname}"
    cube = SpectralCube.read(image)
    position = SkyCoord(ra*u.deg, dec*u.deg)
    size = src_width*u.deg + cube.header['BMAJ']*u.deg*pad
    cutout = cube.subcube(xlo=position.ra - size/2, xhi=position.ra +
                          size/2, ylo=position.dec - size/2, yhi=position.dec + size/2)
    cutout.write(outfile, overwrite=True)
    if verbose:
        print(f'Written to {outfile}')

    image = image.replace('image.restored', 'weights').replace(
        '.total.fits', '.fits')
    outfile = outfile.replace('image.restored', 'weights').replace(
        '.total.fits', '.fits')
    if verbose:
        print(f'Reading {image}')
    cube = SpectralCube.read(image)
    cutout = cube.subcube(xlo=position.ra - size/2, xhi=position.ra +
                          size/2, ylo=position.dec - size/2, yhi=position.dec + size/2)
    cutout.write(outfile, overwrite=True)
    if verbose:
        print(f'Written to {outfile}')

    return src_name


def get_args(args):
    island, beam, island_id, outdir, field, datadir, pad, verbose, dryrun = args

    assert island['Source_ID'] == island_id
    assert beam['Source_ID'] == island_id

    beam_list = np.unique(beam['beams'][field]['beam_list'])

    outdir = f"{outdir}/{island['Source_ID']}"
    try:
        os.mkdir(outdir)
        if verbose:
            print('Made island directory.')
    except FileExistsError:
        if verbose:
            print('Directory exists.')

    args = []
    for beam_num in beam_list:
        images = glob(
            f'{datadir}/image.restored*contcube*beam{beam_num:02}.total.fits')
        for image in images:
            args += [image,
                     island['Source_ID'],
                     island['RA'],
                     island['Dec'],
                     island['Maj']/60/60,
                     outdir,
                     pad,
                     verbose]
    return args


def cutout_islands(field, directory, pool, host, verbose=True, pad=3, verbose_worker=False, dryrun=True):
    if directory[-1] == '/':
        directory = directory[:-1]

    with pymongo.MongoClient(host=host) as client:
        # default connection (ie, local)
        mydb = client['spiceracs']  # Create/open database
        island_col = mydb['islands']  # Create/open collection
        beams_col = mydb['beams']  # Create/open collection

    query = {
        '$and':  [
            {f'beams.{field}': {'$exists': True}},
            {f'beams.{field}.DR1': True}
        ]
    }

    beams = beams_col.find(query).sort('Source_ID')

    island_ids = sorted(beams_col.distinct('Source_ID', query))
    islands = island_col.find(
        {'Source_ID': {'$in': island_ids}}).sort('Source_ID')
    outdir = f'{directory}/cutouts'

    try:
        os.mkdir(outdir)
        print('Made cutout directory.')
    except FileExistsError:
        print('Directory exists.')

    inputs = [[island, beam, island_id, outdir, field, directory, pad,
               verbose_worker, dryrun] for (island_id, island, beam) in zip(island_ids, islands, beams)]

    inputs = inputs[:200]
    args = list(tqdm(pool.imap(get_args, inputs),  # , chunksize=len(inputs)//20),
                     total=len(inputs),
                     disable=(not verbose),
                     desc='Doing cutouts'
                     ))
    embed()
    # flatten list
    # commands = [
    #     item for sublist in commands for subsublist in sublist for item in subsublist
    # ]
    # # if verbose:
    # print(f"I've got {len(commands)} commands to run!")
    # if not dryrun:
    #     script_dir = os.path.dirname(os.path.realpath(__file__))
    #     check_cmd = f"python {script_dir}/check_cutout.py {len(commands)} {outdir}"
    #     print("Run this if you want to follow progress:")
    #     print(check_cmd)
    #     #os.spawnl(os.P_NOWAIT, *shlex.split(check_cmd))
    #     run_command_partial = partial(run_command, verbose=verbose_worker)
    #     failed_commands = list(
    #         tqdm(
    #             pool.imap(run_command_partial, commands),
    #             total=(len(commands)),
    #             disable=(not verbose),
    #             desc='Extracting cubelets'
    #         )
    #     )
    #     real_failed = [command for command in failed_commands if command is not None]
    #     fail_file = f'{directory}/{field}_failedcmds.txt'
    #     if verbose:
    #         print(f"Writing failed cmds to {fail_file}")
    #     with open(fail_file, 'w') as f:
    #         for failed in sorted(real_failed):
    #             f.write(failed+' &'+'\n')
    # TODO: Re-run failed commands
    # It seems to work


def main(args, pool, verbose=True):
    """Main script

    Arguments:
        args {[type]} -- commandline args
    """
    cutout_islands(args.field,
                   args.datadir,
                   pool,
                   args.host,
                   verbose=verbose,
                   pad=args.pad,
                   verbose_worker=args.verbose_worker,
                   dryrun=args.dryrun)

    print("Done!")


def cli():
    """Command-line interface
    """
    import argparse

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
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'datadir',
        metavar='datadir',
        type=str,
        help='Directory containing data cubes in FITS format.')

    parser.add_argument(
        'field',
        metavar='field',
        type=str,
        help='Name of field (e.g. 2132-50A).')

    parser.add_argument(
        'host',
        metavar='host',
        type=str,
        help='Host of mongodb (probably $hostname -i).')

    parser.add_argument(
        "-v",
        dest="verbose",
        action="store_true",
        help="Verbose output [False]."
    )
    parser.add_argument(
        '-p',
        '--pad',
        dest='pad',
        type=float,
        default=3,
        help='Number of beamwidths to pad around source [3].')

    parser.add_argument(
        "-vw",
        dest="verbose_worker",
        action="store_true",
        help="Verbose worker output [False]."
    )
    parser.add_argument(
        "-d",
        dest="dryrun",
        action="store_true",
        help="Do a dry-run [False]."
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
    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)
    if args.mpi:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
    # make it so we can use imap in serial and mpi mode
    if not isinstance(pool, schwimmbad.MultiPool):
        pool.imap = pool.map

    verbose = args.verbose
    if verbose:
        print('Testing MongoDB connection...')
    # default connection (ie, local)
    with pymongo.MongoClient(host=args.host) as client:
        try:
            client.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")
        else:
            if verbose:
                print('MongoDB connection succesful!')

    if verbose:
        print(f"Using pool: {pool.__class__.__name__}")
    main(args, pool, verbose=verbose)
    pool.close()


if __name__ == "__main__":
    cli()
