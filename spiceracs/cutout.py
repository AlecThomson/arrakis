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
from spiceracs.utils import getdata, MyEncoder, try_mkdir, tqdm_dask, get_db, test_db
from astropy.utils import iers
import astropy.units as u
from spectral_cube import SpectralCube
from pprint import pprint

iers.conf.auto_download = False
warnings.filterwarnings(
    "ignore", message="Cube is a Stokes cube, returning spectral cube for I component"
)


@delayed
def cutout(
    image,
    src_name,
    beam,
    ra_hi,
    ra_lo,
    dec_hi,
    dec_lo,
    outdir,
    stoke,
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
    outdir = os.path.abspath(outdir)

    ret = []
    for imtype in ["image", "weight"]:
        basename = os.path.basename(image)
        outname = f"{src_name}.cutout.{basename}"
        outfile = os.path.join(outdir, outname)

        if imtype == "weight":
            image = image.replace("image.restored", "weights").replace(
                ".conv.fits", ".fits"
            )
            outfile = outfile.replace("image.restored", "weights").replace(
                ".conv.fits", ".fits"
            )

        if verbose:
            print(f"Reading {image}")

        cube = SpectralCube.read(image)
        if imtype == "image":
            padder = cube.header["BMAJ"] * u.deg * pad

        xlo = Longitude(ra_lo * u.deg) - Longitude(padder)
        xhi = Longitude(ra_hi * u.deg) + Longitude(padder)
        ylo = Latitude(dec_lo * u.deg) - Latitude(padder)
        yhi = Latitude(dec_hi * u.deg) + Latitude(padder)

        xp_lo, yp_lo = skycoord_to_pixel(SkyCoord(xlo, ylo), cube.wcs)
        xp_hi, yp_hi = skycoord_to_pixel(SkyCoord(xhi, yhi), cube.wcs)

        # Round for cutout
        yp_lo_idx = int(np.floor(yp_lo))
        yp_hi_idx = int(np.ceil(yp_hi))
        xp_lo_idx = int(np.floor(xp_hi))
        xp_hi_idx = int(np.ceil(xp_lo))

        # Use subcube for header transformation
        cutout_cube = cube[:, yp_lo_idx:yp_hi_idx, xp_lo_idx:xp_hi_idx]

        with fits.open(image, memmap=True, mode="denywrite") as hdulist:
            data = hdulist[0].data

            sub_data = data[
                :, 0, yp_lo_idx:yp_hi_idx, xp_lo_idx:xp_hi_idx  # freq  # useless Stokes
            ]
        if not dryrun:
            fits.writeto(outfile, sub_data, header=cutout_cube.header, overwrite=True)
            if verbose:
                print(f"Written to {outfile}")

        # Update database
        myquery = {"Source_ID": src_name}

        filename = os.path.join(
            os.path.basename(os.path.dirname(outfile)), os.path.basename(outfile)
        )
        newvalues = {
            "$set": {f"beams.{field}.{stoke}_beam{beam}_{imtype}_file": filename}
        }

        ret += [pymongo.UpdateOne(myquery, newvalues, upsert=True)]

    return ret


@delayed
def get_args(
    island, comps, beam, island_id, outdir, field, datadir, stokeslist, verbose=True
):
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
    assert island["Source_ID"] == island_id
    assert beam["Source_ID"] == island_id

    beam_list = list(set(beam['beams'][field]['beam_list']))

    outdir = f"{outdir}/{island['Source_ID']}"
    try_mkdir(outdir, verbose=verbose)

    # Find image size
    ras = []
    decs = []
    majs = []
    for comp in comps:
        ras = ras + [comp["RA"]]
        decs = decs + [comp["Dec"]]
        majs = majs + [comp["Maj"]]
    ras = ras * u.deg
    decs = decs * u.deg
    majs = majs * u.arcsec
    coords = SkyCoord(ras, decs)

    try:
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
    except Exception as e:
        print("coords are", coords)
        print("comps are", comps)
        raise e

    args = []
    for beam_num in beam_list:
        for stoke in stokeslist:
            wild = f"{datadir}/image.restored.{stoke.lower()}*contcube*beam{beam_num:02}.conv.fits"
            images = glob(wild)
            if len(images) == 0:
                raise Exception(f"No images found matching '{wild}'")
            for image in images:
                args.extend(
                    [
                        {
                            "image": image,
                            "id": island["Source_ID"],
                            "ra_hi": ra_hi.deg,
                            "ra_lo": ra_lo.deg,
                            "dec_hi": dec_hi.deg,
                            "dec_lo": dec_lo.deg,
                            "outdir": outdir,
                            "beam": beam_num,
                            "stoke": stoke.lower(),
                        }
                    ]
                )
    return args


@delayed
def find_comps(island_id, comp_col):
    comps = list(comp_col.find({"Source_ID": island_id}))
    return comps


@delayed
def unpack(list_sq):
    list_fl = []
    for i in list_sq:
        for j in i:
            list_fl.append(j)
    return list_fl


def cutout_islands(
    field,
    directory,
    host,
    client,
    username=None,
    password=None,
    verbose=True,
    pad=3,
    stokeslist=None,
    verbose_worker=False,
    dryrun=True,
):
    if stokeslist is None:
        stokeslist = ["I", "Q", "U", "V"]
    if verbose:
        print(f"Client is {client}")
    directory = os.path.abspath(directory)
    outdir = os.path.join(directory, "cutouts")

    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )

    # Query the DB
    query = {
        "$and": [{f"beams.{field}": {"$exists": True}}, {f"beams.{field}.DR1": True}]
    }

    beams = list(beams_col.find(query).sort("Source_ID"))

    island_ids = sorted(beams_col.distinct("Source_ID", query))
    islands = list(island_col.find({"Source_ID": {"$in": island_ids}}).sort("Source_ID"))

    big_comps = list(comp_col.find({"Source_ID": {'$in': island_ids}}))
    comps = []
    for island_id in island_ids:
        _comps = []
        for c in big_comps:
            if c["Source_ID"] == island_id:
                _comps.append(c)
        comps.append(_comps)

    # Create output dir if it doesn't exist
    try_mkdir(outdir)

    args = []
    for (island_id, island, comp, beam) in zip(island_ids, islands, comps, beams):
        if len(comp) == 0:
            warnings.warn(f"Skipping island {island_id} -- no components found")
            continue
        else:
            arg = get_args(
                island,
                comp,
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
    tqdm_dask(
        flat_args, desc="Getting args", disable=(not verbose), total=len(islands)+1
    )
    flat_args = flat_args.result()
    cuts = []
    for arg in flat_args:
        cut = cutout(
            image=arg["image"],
            src_name=arg["id"],
            beam=arg["beam"],
            ra_hi=arg["ra_hi"],
            ra_lo=arg["ra_lo"],
            dec_hi=arg["dec_hi"],
            dec_lo=arg["dec_lo"],
            outdir=arg["outdir"],
            stoke=arg["stoke"],
            field=field,
            pad=pad,
            verbose=verbose_worker,
            dryrun=dryrun,
        )
        cuts.append(cut)

    futures = client.persist(cuts)
    # dumb solution for https://github.com/dask/distributed/issues/4831
    time.sleep(5)
    tqdm_dask(futures, desc="Cutting out", disable=(not verbose))
    if not dryrun:
        _updates = [f.compute() for f in futures]
        updates = [val for sublist in _updates for val in sublist]
        if verbose:
            print("Updating database...")
        db_res = beams_col.bulk_write(updates, ordered=False)
        if verbose:
            pprint(db_res.bulk_api_result)

    print("Cutouts Done!")

def main(args, verbose=True):
    """Main script

    Arguments:
        args {[type]} -- commandline args
    """
    cluster = LocalCluster(
        n_workers=12, threads_per_worker=1, dashboard_address=":9898"
    )
    client = Client(cluster)
    print(client)
    cutout_islands(
        field=args.field,
        directory=args.datadir,
        host=args.host,
        client=client,
        username=args.username,
        password=args.password,
        verbose=verbose,
        pad=args.pad,
        stokeslist=args.stokeslist,
        verbose_worker=args.verbose_worker,
        dryrun=args.dryrun,
    )

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
        "field", metavar="field", type=str, help="Name of field (e.g. 2132-50A)."
    )

    parser.add_argument(
        "datadir",
        metavar="datadir",
        type=str,
        help="Directory containing data cubes in FITS format.",
    )

    parser.add_argument(
        "host",
        metavar="host",
        type=str,
        help="Host of mongodb (probably $hostname -i).",
    )

    parser.add_argument(
        "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose output [False]."
    )
    parser.add_argument(
        "-p",
        "--pad",
        dest="pad",
        type=float,
        default=3,
        help="Number of beamwidths to pad around source [3].",
    )

    parser.add_argument(
        "-vw",
        dest="verbose_worker",
        action="store_true",
        help="Verbose worker output [False].",
    )
    parser.add_argument(
        "-d", "--dryrun", action="store_true", help="Do a dry-run [False]."
    )
    parser.add_argument(
        "-s",
        "--stokes",
        dest="stokeslist",
        nargs="+",
        type=str,
        help="List of Stokes parameters to image [ALL]",
    )

    args = parser.parse_args()

    verbose = args.verbose
    test_db(
        host=args.host,
        username=args.username,
        password=args.password,
        verbose=verbose
    )

    main(args, verbose=verbose)


if __name__ == "__main__":
    cli()
