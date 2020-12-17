#!/usr/bin/env python
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
iers.conf.auto_download = False
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print = functools.partial(print, f'[{psutil.Process().cpu_num()},{rank}]', flush=True)


def cut_command(image, src_name, ra, dec, src_width, outdir, pad=3, verbose=False):
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
    if outdir[-1] == '/':
        outdir = outdir[:-1]

    with fits.open(image, memmap=True, mode='denywrite') as hdulist:
        hdu = hdulist[0]
        header = hdu.header
        shape = hdu.data.shape
    wcs = WCS(header)
    freq = 888e6
    x_s, y_s, _, _ = np.array(wcs.all_world2pix(
        ra, dec, 0, freq, 0)).astype(int)

    pixels_per_beam = int(header['BMAJ']/header['CDELT2'])

    src_width_pix = int(src_width / header['CDELT2'])

    y_max = y_s + src_width_pix
    x_max = x_s + src_width_pix
    y_min = y_s - src_width_pix
    x_min = x_s - src_width_pix

    # Skip if source is outside of cube bounds
    if (y_max > header['NAXIS2'] or x_max > header['NAXIS1'] or
            x_min < 0 or y_min < 0):
        return

    starty = int(y_min-pad*pixels_per_beam)
    stopy = int(y_max+pad*pixels_per_beam)
    startx = int(x_min-pad*pixels_per_beam)
    stopx = int(x_max+pad*pixels_per_beam)

    # Check if pad puts bbox outside of cube
    if starty < 0:
        starty = y_min
    if startx < 0:
        startx = x_min
    if stopy > header['NAXIS2']:
        stopy = y_max
    if stopx > header['NAXIS1']:
        stopx = x_max

    basename = os.path.basename(image)
    outname = f'{src_name}.cutout.{basename}'
    outfile = f"{outdir}/{outname}"

    #command_dict = {}

    command_image = f"fitscopy '{image}[{startx+1}:{stopx},{starty+1}:{stopy}]' '!{outfile}'"

    # if verbose:
    #    print(command)
    #command = shlex.split(command)

    # command_dict.update(
    #    {
    #        'image': command
    #    }
    # )

    # Now do weights
    image = image.replace('image.restored', 'weights').replace('.total.fits','.fits')
    outfile = outfile.replace('image.restored', 'weights').replace('.total.fits','.fits')

    command_weight = f"fitscopy '{image}[{startx+1}:{stopx},{starty+1}:{stopy}]' '!{outfile}'"

    #command = shlex.split(command)
    #
    # command_dict.update(
    #    {
    #        'weight': command
    #    }
    # )
    return [command_image, command_weight]


def run_command(command, verbose=False):
    try:
        proc = subprocess.run(shlex.split(command),
                              stderr=subprocess.STDOUT,
                              encoding='utf-8'
                              )

        retries = 0
        while proc.returncode != 0:
            proc = subprocess.run(shlex.split(command),
                                  stderr=subprocess.STDOUT,
                                  encoding='utf-8'
                                  )
            retries += 1

            if retries > 1e6:
                break

        if proc.returncode == 0:
            if verbose:
                print(command)
            return

        elif proc.returncode != 0:
            return command

    except:
        print('I failed in my job!', command)


def get_cut_command(args):
    island, beam, island_id, outdir, field, datadir, pad, verbose, dryrun = args

    assert island['Source_ID'] == island_id
    assert beam['Source_ID'] == island_id

    beam_list = beam['beams'][field]['beam_list']

    outdir = f"{outdir}/{island['Source_ID']}"
    try:
        os.mkdir(outdir)
        if verbose:
            print('Made island directory.')
    except FileExistsError:
        if verbose:
            print('Directory exists.')

    commands = []
    for beam_num in beam_list:
        images = glob(
            f'{datadir}/image.restored*contcube*beam{beam_num:02}.total.fits')
        for image in images:
            command = cut_command(image,
                                  island['Source_ID'],
                                  island['RA'],
                                  island['Dec'],
                                  island['Maj']/60/60,
                                  outdir,
                                  pad=pad, verbose=verbose)
            commands.append(command)
    return commands


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
    islands = island_col.find({'Source_ID': {'$in': island_ids}}).sort('Source_ID')
    outdir = f'{directory}/cutouts'

    try:
        os.mkdir(outdir)
        print('Made cutout directory.')
    except FileExistsError:
        print('Directory exists.')

    inputs = [[island, beam, island_id, outdir, field, directory, pad,
               verbose_worker, dryrun] for (island_id, island, beam) in zip(island_ids, islands, beams)]

    commands = list(tqdm(pool.imap(get_cut_command, inputs),  # , chunksize=len(inputs)//20),
                         total=len(inputs),
                         disable=(not verbose),
                         desc='Generating commands'
                         ))
    # flatten list
    commands = [
        item for sublist in commands for subsublist in sublist for item in subsublist
    ]
    # if verbose:
    print(f"I've got {len(commands)} commands to run!")
    if not dryrun:
        script_dir = os.path.dirname(os.path.realpath(__file__))
        check_cmd = f"python {script_dir}/check_cutout.py {len(commands)} {outdir}"
        print("Run this if you want to follow progress:")
        print(check_cmd)
        #os.spawnl(os.P_NOWAIT, *shlex.split(check_cmd))
        run_command_partial = partial(run_command, verbose=verbose_worker)
        failed_commands = list(
            tqdm(
                pool.imap(run_command_partial, commands),
                total=(len(commands)),
                disable=(not verbose),
                desc='Extracting cubelets'
            )
        )
        real_failed = [command for command in failed_commands if command is not None]
        fail_file = f'{directory}/{field}_failedcmds.txt'
        if verbose:
            print(f"Writing failed cmds to {fail_file}")
        with open(fail_file, 'w') as f: 
            for failed in sorted(real_failed):
                f.write(failed+' &'+'\n') 
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
    client = pymongo.MongoClient(host=args.host)
    try:
        client.list_database_names()
    except pymongo.errors.ServerSelectionTimeoutError:
        raise Exception("Please ensure 'mongod' is running")
    else:
        if verbose:
            print('MongoDB connection succesful!')
    client.close()

    if verbose:
        print(f"Using pool: {pool.__class__.__name__}")
    main(args, pool, verbose=verbose)
    pool.close()


if __name__ == "__main__":
    cli()
