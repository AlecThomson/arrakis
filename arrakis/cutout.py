#!/usr/bin/env python
"""Produce cutouts from RACS cubes"""
import argparse
import logging
import os
import warnings
from glob import glob
from pprint import pformat
from shutil import copyfile
from typing import Dict, List, Union

import astropy.units as u
import numpy as np
import pymongo
from astropy import units as u
from astropy.coordinates import Latitude, Longitude, SkyCoord
from astropy.io import fits
from astropy.utils import iers
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs.utils import skycoord_to_pixel
from dask import delayed
from dask.distributed import Client, LocalCluster, progress
from distributed import get_client
from IPython import embed
from spectral_cube import SpectralCube
from spectral_cube.utils import SpectralCubeWarning

from arrakis.logger import logger
from arrakis.utils.database import get_db, test_db
from arrakis.utils.fitsutils import fix_header
from arrakis.utils.io import try_mkdir
from arrakis.utils.pipeline import chunk_dask, logo_str, tqdm_dask

iers.conf.auto_download = False
warnings.filterwarnings(
    "ignore", message="Cube is a Stokes cube, returning spectral cube for I component"
)

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)
warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")

logger.setLevel(logging.INFO)


@delayed
def cutout(
    image: str,
    src_name: str,
    beam: int,
    ra_hi: float,
    ra_lo: float,
    dec_hi: float,
    dec_lo: float,
    outdir: str,
    stoke: str,
    field: str,
    pad=3,
    verbose=False,
    dryrun=False,
) -> List[pymongo.UpdateOne]:
    """Perform a cutout.

    Args:
        image (str): Name of the image file
        src_name (str): Name of the RACS source
        beam (int): Beam number
        ra_hi (float): Upper RA bound
        ra_lo (float): Lower RA bound
        dec_hi (float): Upper DEC bound
        dec_lo (float): Lower DEC bound
        outdir (str): Output directgory
        stoke (str): Stokes parameter
        field (str): RACS field name
        pad (int, optional): Number of beamwidths to pad. Defaults to 3.
        verbose (bool, optional): Verbose output. Defaults to False.
        dryrun (bool, optional): Don't save FITS files. Defaults to False.

    Returns:
        pymongo.UpdateOne: Update query for MongoDB
    """
    logger.setLevel(logging.INFO)
    # logger = logging.getLogger('distributed.worker')
    # logger = get_run_logger()

    logger.info(f"Timwashere - {image=}")

    outdir = os.path.abspath(outdir)

    ret = []
    for imtype in ["image", "weight"]:
        basename = os.path.basename(image)
        outname = f"{src_name}.cutout.{basename}"
        outfile = os.path.join(outdir, outname)

        if imtype == "weight":
            image = image.replace("image.restored", "weights.restored").replace(
                ".fits", ".txt"
            )
            outfile = outfile.replace("image.restored", "weights.restored").replace(
                ".fits", ".txt"
            )
            copyfile(image, outfile)
            logger.info(f"Written to {outfile}")

        if imtype == "image":
            logger.info(f"Reading {image}")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", AstropyWarning)
                cube = SpectralCube.read(image)
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
            new_header = cutout_cube.header
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", AstropyWarning)
                with fits.open(image, memmap=True, mode="denywrite") as hdulist:
                    data = hdulist[0].data
                    old_header = hdulist[0].header

                    sub_data = data[
                        :,
                        :,
                        yp_lo_idx:yp_hi_idx,
                        xp_lo_idx:xp_hi_idx,  # freq, Stokes, y, x
                    ]
                fixed_header = fix_header(new_header, old_header)
                # Add source name to header for CASDA
                fixed_header["OBJECT"] = src_name
                if not dryrun:
                    fits.writeto(
                        outfile,
                        sub_data,
                        header=fixed_header,
                        overwrite=True,
                        output_verify="fix",
                    )
                    logger.info(f"Written to {outfile}")

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
    island: Dict,
    comps: List[Dict],
    beam: Dict,
    island_id: str,
    outdir: str,
    field: str,
    datadir: str,
    stokeslist: List[str],
    verbose=True,
) -> List[Dict]:
    """Get arguments for cutout function

    Args:
        island (str): Mongo entry for RACS island
        comps (List[Dict]): List of mongo entries for RACS components in island
        beam (Dict): Mongo entry for the RACS beam
        island_id (str): RACS island ID
        outdir (str): Output directory
        field (str): RACS field name
        datadir (str): Input directory
        stokeslist (List[str]): List of Stokes parameters to process
        verbose (bool, optional): Verbose output. Defaults to True.

    Raises:
        e: Exception
        Exception: Problems with coordinates

    Returns:
        List[Dict]: List of cutout arguments for cutout function
    """

    logger.setLevel(logging.INFO)

    assert island["Source_ID"] == island_id
    assert beam["Source_ID"] == island_id

    beam_list = list(set(beam["beams"][field]["beam_list"]))

    outdir = f"{outdir}/{island['Source_ID']}"
    try_mkdir(outdir, verbose=verbose)

    # Find image size
    ras = []  # type: List[float]
    decs = []  # type: List[float]
    majs = []  # type: List[float]
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
        logger.debug(f"coords are {coords=}")
        logger.debug(f"comps are {comps=}")
        raise e

    args = []
    for beam_num in beam_list:
        for stoke in stokeslist:
            wild = f"{datadir}/image.restored.{stoke.lower()}*contcube*beam{beam_num:02}.conv.fits"
            images = glob(wild)
            if len(images) == 0:
                raise Exception(f"No images found matching '{wild}'")
            elif len(images) > 1:
                raise Exception(
                    f"More than one image found matching '{wild}'. Files {images=}"
                )

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
def find_comps(island_id: str, comp_col: pymongo.collection.Collection) -> List[Dict]:
    """Find components for a given island

    Args:
        island_id (str): RACS island ID
        comp_col (pymongo.collection.Collection): Component collection

    Returns:
        List[Dict]: List of mongo entries for RACS components in island
    """
    comps = list(comp_col.find({"Source_ID": island_id}))
    return comps


@delayed
def unpack(list_sq: List[List[Dict]]) -> List[Dict]:
    """Unpack list of lists

    Args:
        list_sq (List[List[Dict]]): List of lists of dicts

    Returns:
        List[Dict]: List of dicts
    """
    list_fl = []
    for i in list_sq:
        for j in i:
            list_fl.append(j)
    return list_fl


def cutout_islands(
    field: str,
    directory: str,
    host: str,
    username: Union[str, None] = None,
    password: Union[str, None] = None,
    pad: float = 3,
    stokeslist: Union[List[str], None] = None,
    verbose_worker: bool = False,
    dryrun: bool = True,
) -> None:
    """Perform cutouts of RACS islands in parallel.

    Args:
        field (str): RACS field name.
        directory (str): Directory to store cutouts.
        host (str): MongoDB host.
        client (Client): Dask client.
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo password. Defaults to None.
        verbose (bool, optional): Verbose output. Defaults to True.
        pad (int, optional): Number of beamwidths to pad cutouts. Defaults to 3.
        stokeslist (List[str], optional): Stokes parameters to cutout. Defaults to None.
        verbose_worker (bool, optional): Worker function outout. Defaults to False.
        dryrun (bool, optional): Do everything except write FITS files. Defaults to True.
    """
    if stokeslist is None:
        stokeslist = ["I", "Q", "U", "V"]
    client = get_client()
    logger.debug(f"Client is {client}")
    directory = os.path.abspath(directory)
    outdir = os.path.join(directory, "cutouts")

    logger.info("Testing database. ")
    test_db(
        host=host,
        username=username,
        password=password,
    )

    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )

    # Query the DB
    query = {"$and": [{f"beams.{field}": {"$exists": True}}]}

    beams = list(beams_col.find(query).sort("Source_ID"))

    island_ids = sorted(beams_col.distinct("Source_ID", query))
    islands = list(
        island_col.find({"Source_ID": {"$in": island_ids}}).sort("Source_ID")
    )

    big_comps = list(comp_col.find({"Source_ID": {"$in": island_ids}}))
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
    for island_id, island, comp, beam in zip(island_ids, islands, comps, beams):
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
    tqdm_dask(flat_args, desc="Getting args", total=len(islands) + 1)
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

    futures = chunk_dask(
        outputs=cuts,
        task_name="cutouts",
        progress_text="Cutting out",
    )
    if not dryrun:
        _updates = [f.compute() for f in futures]
        updates = [val for sublist in _updates for val in sublist]
        logger.info("Updating database...")
        db_res = beams_col.bulk_write(updates, ordered=False)
        logger.info(pformat(db_res.bulk_api_result))

    logger.info("Cutouts Done!")


def main(args: argparse.Namespace, verbose=True) -> None:
    """Main script

    Args:
        args (argparse.Namespace): Command-line args
        verbose (bool, optional): Verbose output. Defaults to True.
    """
    cutout_islands(
        field=args.field,
        directory=args.datadir,
        host=args.host,
        username=args.username,
        password=args.password,
        pad=args.pad,
        stokeslist=args.stokeslist,
        verbose_worker=args.verbose_worker,
        dryrun=args.dryrun,
    )

    logger.info("Done!")


def cutout_parser(parent_parser: bool = False) -> argparse.ArgumentParser:
    descStr = f"""
    {logo_str}

    Arrakis Stage 1:
    Produce cubelets from a RACS field using a Selavy table.
    If Stokes V is present, it will be squished into RMS spectra.

    To use with MPI:
       mpirun -n $NPROCS python -u cutout.py $cubedir $tabledir
       $outdir --mpi
    """

    # Parse the command line options
    cut_parser = argparse.ArgumentParser(
        add_help=not parent_parser,
        description=descStr,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser = cut_parser.add_argument_group("cutout arguments")

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

    return cut_parser


def cli() -> None:
    """Command-line interface"""
    parser = cutout_parser()

    args = parser.parse_args()

    verbose = args.verbose
    if verbose:
        logger.setLevel(logging.INFO)

    cluster = LocalCluster(
        n_workers=12, threads_per_worker=1, dashboard_address=":9898"
    )
    client = Client(cluster)
    logger.info(client)

    test_db(
        host=args.host,
        username=args.username,
        password=args.password,
    )

    main(args, verbose=verbose)


if __name__ == "__main__":
    cli()
