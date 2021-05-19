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
import json
import time
from tqdm import tqdm, trange
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import units as u
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy.io import fits
import sys
import os
from glob import glob
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord, Longitude, Latitude
from spiceracs.utils import getdata, MyEncoder
from astropy.utils import iers
import astropy.units as u
from spectral_cube import SpectralCube
iers.conf.auto_download = False
warnings.filterwarnings(
    "ignore", message="Cube is a Stokes cube, returning spectral cube for I component")


@delayed
def cutout(image,
           src_name,
           beam,
           ra_hi,
           ra_lo,
           dec_hi,
           dec_lo,
           outdir,
           stoke,
           host,
           field,
           pad=3,
           verbose=False,
           dryrun=False
           ):
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
    with pymongo.MongoClient(host=host, connect=False) as dbclient:
        mydb = dbclient['spiceracs']  # Create/open database
        beams_col = mydb['beams']  # Create/open collection

    if verbose:
        print(f'Reading {image}')
    if outdir[-1] == '/':
        outdir = outdir[:-1]
    basename = os.path.basename(image)
    outname = f'{src_name}.cutout.{basename}'
    outfile = f"{outdir}/{outname}"
    cube = SpectralCube.read(image)
    padder = cube.header['BMAJ']*u.deg * pad

    # Don't know why this is required
    # Dask breaks without it...
    # ra_lo = ra_lo.value*u.arcsec
    # ra_hi = ra_hi.value*u.arcsec
    # dec_lo = dec_lo.value*u.arcsec
    # dec_hi = dec_hi.value*u.arcsec

    # xlo = np.sign(ra_lo) * (np.abs(ra_lo) + padder) #ra_lo - padder
    # xhi = np.sign(ra_hi) * (np.abs(ra_hi) + padder) #ra_hi + padder
    # ylo = np.sign(dec_lo) * (np.abs(dec_lo) + padder) #dec_lo  - padder
    # yhi = np.sign(dec_hi) * (np.abs(dec_hi) + padder) #dec_hi  + padder

    xlo = Longitude(ra_lo*u.deg) - Longitude(padder)
    xhi = Longitude(ra_hi*u.deg) + Longitude(padder)
    ylo = Latitude(dec_lo*u.deg) - Latitude(padder)
    yhi = Latitude(dec_hi*u.deg) + Latitude(padder)

    xp_lo, yp_lo = skycoord_to_pixel(SkyCoord(xlo, ylo), cube.wcs)
    xp_hi, yp_hi = skycoord_to_pixel(SkyCoord(xhi, yhi), cube.wcs)

    cutout_cube = cube[
        :,
        int(np.floor(yp_lo)):int(np.ceil(yp_hi)),
        int(np.floor(xp_hi)):int(np.ceil(xp_lo))
    ]

    # cutout_cube = cube.subcube(xlo=xlo.deg*u.deg,
    #                       xhi=xhi.deg*u.deg,
    #                       ylo=ylo.deg*u.deg,
    #                       yhi=yhi.deg*u.deg,
    #                       )
    if not dryrun:
        fits.writeto(outfile,
                     cutout_cube.unmasked_data[:].value,
                     header=cutout_cube.header,
                     overwrite=True
                     )
        # cutout_cube.write(outfile, overwrite=True)
        if verbose:
            print(f'Written to {outfile}')

        # Update database
        myquery = {
            "Source_ID": src_name
        }
        newvalues = {
            "$set":
            {
                f"beams.{field}.{stoke}_beam{beam}_image_file": outfile
            }
        }

        beams_col.update_one(myquery, newvalues, upsert=True)

    image = image.replace('image.restored', 'weights').replace(
        '.conv.fits', '.fits')
    outfile = outfile.replace('image.restored', 'weights').replace(
        '.conv.fits', '.fits')
    if verbose:
        print(f'Reading {image}')
    cube = SpectralCube.read(image)
    cutout_cube = cube[
        :,
        int(np.floor(yp_lo)):int(np.ceil(yp_hi)),
        int(np.floor(xp_hi)):int(np.ceil(xp_lo))
    ]
    # cutout_cube = cube.subcube(xlo=xlo.deg*u.deg,
    #                       xhi=xhi.deg*u.deg,
    #                       ylo=ylo.deg*u.deg,
    #                       yhi=yhi.deg*u.deg,
    #                       )
    if not dryrun:
        fits.writeto(outfile,
                     cutout_cube.unmasked_data[:].value,
                     header=cutout_cube.header,
                     overwrite=True
                     )
        # cutout_cube.write(outfile, overwrite=True)
        if verbose:
            print(f'Written to {outfile}')

        # Update database
        myquery = {
            "Source_ID": src_name
        }
        newvalues = {
            "$set":
            {
                f"beams.{field}.{stoke}_beam{beam}_weight_file": outfile
            }
        }

        beams_col.update_one(myquery, newvalues, upsert=True)

    return src_name


@delayed
def get_args(island,
             comps,
             beam,
             island_id,
             outdir,
             field,
             datadir,
             stokeslist,
             verbose=True):
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
        ras = ras + [comp['RA']]
        decs = decs + [comp['Dec']]
        majs = majs + [comp['Maj']]
    ras = ras * u.deg
    decs = decs * u.deg
    majs = majs * u.arcsec
    coords = SkyCoord(ras, decs)

    # ra_hi = (majs[np.argmax(ras)] + np.max(ras)).to(u.arcsec)
    # ra_lo = (np.min(ras) - majs[np.argmin(ras)]).to(u.arcsec)
    # dec_hi = (majs[np.argmax(decs)] + np.max(decs)).to(u.arcsec)
    # dec_lo = (np.min(decs) - majs[np.argmin(decs)]).to(u.arcsec)

    ra_max = np.max(coords.ra)
    ra_i_max = np.argmax(coords.ra)
    ra_off = Longitude(majs[ra_i_max])
    ra_hi = ra_max + ra_off

    ra_min = np.min(coords.ra)
    ra_i_min = np.argmin(coords.ra)
    ra_off = Longitude(majs[ra_i_min])
    ra_lo = ra_min - ra_off

    dec_max = np.max(coords.dec)
    dec_i_max = np.argmax(coords.dec)
    dec_off = Longitude(majs[dec_i_max])
    dec_hi = dec_max + dec_off

    dec_min = np.min(coords.dec)
    dec_i_min = np.argmin(coords.dec)
    dec_off = Longitude(majs[dec_i_min])
    dec_lo = dec_min - dec_off

    args = []
    for beam_num in beam_list:
        for stoke in stokeslist:
            wild = f'{datadir}/image.restored.{stoke.lower()}*contcube*beam{beam_num:02}.conv.fits'
            images = glob(
                wild
            )
            if len(images) == 0:
                raise Exception(
                    f"No images found matching '{wild}'"
                )
            for image in images:
                args.extend([
                    {
                        'image': image,
                        'id': island['Source_ID'],
                        'ra_hi': ra_hi.deg,
                        'ra_lo': ra_lo.deg,
                        'dec_hi': dec_hi.deg,
                        'dec_lo': dec_lo.deg,
                        'outdir': outdir,
                        'beam': beam_num,
                        'stoke': stoke.lower()
                    }
                ]
                )
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


def cutout_islands(field,
                   directory,
                   host,
                   client,
                   verbose=True,
                   pad=3,
                   stokeslist=None,
                   verbose_worker=False,
                   dryrun=True):
    if stokeslist is None:
        stokeslist = ["I", "Q", "U", "V"]
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
            stokeslist,
            verbose=verbose_worker,
        )
        args.append(arg)

    flat_args = unpack(args)
    flat_args = client.compute(flat_args)
    if verbose:
        print("Getting args...")
    progress(flat_args)
    flat_args = flat_args.result()
    cuts = []
    for arg in flat_args:
        cut = cutout(
            image=arg['image'],
            src_name=arg['id'],
            beam=arg['beam'],
            ra_hi=arg['ra_hi'],
            ra_lo=arg['ra_lo'],
            dec_hi=arg['dec_hi'],
            dec_lo=arg['dec_lo'],
            outdir=arg['outdir'],
            stoke=arg['stoke'],
            host=host,
            field=field,
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
    cluster = LocalCluster(n_workers=20, dashboard_address=":9999")
    client = Client(cluster)
    cutout_islands(args.field,
                   args.datadir,
                   args.host,
                   client,
                   verbose=verbose,
                   pad=args.pad,
                   stokeslist=args.stokeslist,
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
    parser.add_argument(
        "-s",
        "--stokes",
        dest="stokeslist",
        nargs='+',
        type=str,
        help="List of Stokes parameters to image [ALL]",
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
