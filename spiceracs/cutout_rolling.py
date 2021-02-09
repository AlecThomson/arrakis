#!/usr/bin/env python
import warnings
import pymongo
import dask
from dask import delayed
from dask.distributed import Client, progress, LocalCluster
from dask.diagnostics import ProgressBar
from IPython import embed
from functools import partial
import functools
import psutil
import shlex
import subprocess
import schwimmbad
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
import astropy.units as u
from spectral_cube import SpectralCube
iers.conf.auto_download = False
warnings.filterwarnings(
    "ignore", message="Cube is a Stokes cube, returning spectral cube for I component")


@delayed
def cutout(image, src_name, ra_hi, ra_lo, dec_hi, dec_lo, outdir, pad=3, verbose=False, dryrun=False):
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
    padder = cube.header['BMAJ']*u.deg * pad
    cutout = cube.subcube(xlo=ra_lo - padder,
                          xhi=ra_hi + padder,
                          ylo=dec_lo - padder,
                          yhi=dec_hi + padder
                          )
    if not dryrun:
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
    cutout = cube.subcube(xlo=ra_lo - padder,
                          xhi=ra_hi + padder,
                          ylo=dec_lo - padder,
                          yhi=dec_hi + padder
                          )
    if not dryrun:
        cutout.write(outfile, overwrite=True)
        if verbose:
            print(f'Written to {outfile}')

    return src_name


@delayed
def get_args(island, comps, beam, island_id, outdir, field, datadir, verbose=True):
    """[summary]

    Args:
        island (bool): [description]
        beam ([type]): [description]
        island_id (bool): [description]
        outdir ([type]): [description]
        field ([type]): [description]
        datadir ([type]): [description]
        verbose ([type]): [description]
        dryrun ([type]): [description]

    Returns:
        [type]: [description]
    """
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

    # Find image size
    ras = []
    decs = []
    majs = []
    for comp in comps:
        ras += [comp['RA']]
        decs += [comp['Dec']]
        majs += [comp['Maj']]
    ras *= u.deg
    decs *= u.deg
    majs *= u.arcsec
    ra_hi = majs[np.argmax(ras)] + np.max(ras)
    ra_lo = majs[np.argmin(ras)] + np.min(ras)
    dec_hi = majs[np.argmax(decs)] + np.max(decs)
    dec_lo = majs[np.argmin(decs)] + np.min(decs)

    args = []
    for beam_num in beam_list:
        images = glob(
            f'{datadir}/image.restored*contcube*beam{beam_num:02}.total.fits')
        for image in images:
            args += [
                {
                    'image': image,
                    'id': island['Source_ID'],
                    'ra_hi': ra_hi,
                    'ra_lo': ra_lo,
                    'dec_hi': dec_hi,
                    'dec_lo': dec_lo,
                    'outdir': outdir,
                }
            ]
    return args


@delayed
def find_comps(island_id, host):
    with pymongo.MongoClient(host=host, connect=False) as dbclient:
        # default connection (ie, local)
        mydb = dbclient['spiceracs']  # Create/open database
        comp_col = mydb['components']  # Create/open collection
    comps = comp_col.find({'Source_ID': island_id})
    comps = [c for c in comps]
    return comps


@delayed
def unpack(list_sq):
    list_fl = []
    for i in list_sq:
        for j in i:
            list_fl.append(j)
    return list_fl


# @delayed
# def cutout_loop(flat_args, pad=5, dryrun=False, verbose=True):
#     cuts = []
#     for arg in tqdm(flat_args,
#                     desc=f'Cutting out',
#                     disable=(not verbose)
#                     ):
#         #for a in arg:
#         cut = cutout(
#             arg['image'],
#             arg['id'],
#             arg['ra_hi'],
#             arg['ra_lo'],
#             arg['dec_hi'],
#             arg['dec_lo'],
#             arg['outdir'],
#             pad=pad,
#             verbose=verbose,
#             dryrun=dryrun
#         )
#         cuts.append(cut)
#     return cuts

def cutout_islands(field, directory, host, client, verbose=True, pad=3, verbose_worker=False, dryrun=True):
    if verbose:
        print(f"Client is {client}")
    if directory[-1] == '/':
        directory = directory[:-1]
    outdir = f'{directory}/cutouts'
    with pymongo.MongoClient(host=host, connect=False) as dbclient:
        # default connection (ie, local)
        mydb = dbclient['spiceracs']  # Create/open database
        comp_col = mydb['components']  # Create/open collection
        island_col = mydb['islands']  # Create/open collection
        beams_col = mydb['beams']  # Create/open collection

    # Query the DB
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

    beams = [i for i in beams]
    islands = [i for i in islands]

    # Create output dir if it doesn't exist
    try:
        os.mkdir(outdir)
        print('Made cutout directory.')
    except FileExistsError:
        print('Directory exists.')

    args = []
    for (island_id, island, beam) in zip(island_ids, islands, beams):
        #comps = comp_col.find({'Source_ID': island_id})
        comps = find_comps(island_id, host)
        arg = get_args(
            island,
            comps,
            beam,
            island_id,
            outdir,
            field,
            directory,
            verbose=verbose_worker,
        )
        args.append(arg)

    flat_args = unpack(args)
    flat_args = flat_args.persist()
    if verbose:
        print("Getting args...")
    progress(flat_args)
    flat_args = flat_args.compute()

    # cuts = cutout_loop(flat_args, pad=pad, dryrun=dryrun, verbose=verbose_worker)
    cuts = []
    for arg in flat_args:
        # for a in arg:
        cut = cutout(
            arg['image'],
            arg['id'],
            arg['ra_hi'],
            arg['ra_lo'],
            arg['dec_hi'],
            arg['dec_lo'],
            arg['outdir'],
            pad=pad,
            verbose=verbose_worker,
            dryrun=dryrun
        )
        cuts.append(cut)

    output = client.persist(cuts)
    if verbose:
        print("Cutting out...")
    progress(output)


def main(args, verbose=True):
    """Main script

    Arguments:
        args {[type]} -- commandline args
    """
    cluster = LocalCluster(n_workers=20)
    client = Client(cluster)
    cutout_islands(args.field,
                   args.datadir,
                   args.host,
                   client,
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
        'field',
        metavar='field',
        type=str,
        help='Name of field (e.g. 2132-50A).')


    parser.add_argument(
        'datadir',
        metavar='datadir',
        type=str,
        help='Directory containing data cubes in FITS format.')

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

    verbose = args.verbose
    if verbose:
        print('Testing MongoDB connection...')
    # default connection (ie, local)
    with pymongo.MongoClient(host=args.host, connect=False) as client:
        try:
            client.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")
        else:
            if verbose:
                print('MongoDB connection succesful!')

    main(args, verbose=verbose)


if __name__ == "__main__":
    cli()
