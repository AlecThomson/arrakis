#!/usr/bin/env python
"""Produce cutouts from RACS cubes"""
import argparse
import logging
import os
import pickle
import warnings
from concurrent.futures import ThreadPoolExecutor
from glob import glob
from pprint import pformat
from shutil import copyfile
from typing import Dict, List
from typing import NamedTuple as Struct
from typing import Optional, Set, TypeVar, Union

import astropy.units as u
import numpy as np
import pymongo
from astropy.coordinates import Latitude, Longitude, SkyCoord
from astropy.io import fits
from astropy.utils import iers
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs.utils import skycoord_to_pixel
from prefect import flow, task, unmapped
from spectral_cube import SpectralCube
from spectral_cube.utils import SpectralCubeWarning
from tqdm.auto import tqdm

from arrakis.logger import TqdmToLogger, UltimateHelpFormatter, logger
from arrakis.utils.database import get_db, test_db
from arrakis.utils.fitsutils import fix_header
from arrakis.utils.io import try_mkdir
from arrakis.utils.pipeline import logo_str

iers.conf.auto_download = False
warnings.filterwarnings(
    "ignore", message="Cube is a Stokes cube, returning spectral cube for I component"
)

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)
warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")

logger.setLevel(logging.INFO)
TQDM_OUT = TqdmToLogger(logger, level=logging.INFO)

T = TypeVar("T")


class CutoutArgs(Struct):
    """Arguments for cutout function"""

    """Name of the source"""
    ra_high: float
    """Upper RA bound in degrees"""
    ra_low: float
    """Lower RA bound in degrees"""
    dec_high: float
    """Upper DEC bound in degrees"""
    dec_low: float
    """Lower DEC bound in degrees"""
    outdir: str
    """Output directory"""


def cutout_weight(
    image_name: str,
    source_id: str,
    cutout_args: CutoutArgs,
    field: str,
    stoke: str,
    beam_num: int,
    dryrun=False,
) -> pymongo.UpdateOne:
    outdir = os.path.abspath(cutout_args.outdir)
    basename = os.path.basename(image_name)
    outname = f"{source_id}.cutout.{basename}"
    outfile = os.path.join(outdir, outname)
    image = image_name.replace("image.restored", "weights.restored").replace(
        ".fits", ".txt"
    )
    outfile = outfile.replace("image.restored", "weights.restored").replace(
        ".fits", ".txt"
    )

    if not dryrun:
        copyfile(image, outfile)
        logger.info(f"Written to {outfile}")

    # Update database
    myquery = {"Source_ID": source_id}

    filename = os.path.join(
        os.path.basename(os.path.dirname(outfile)), os.path.basename(outfile)
    )
    newvalues = {
        "$set": {f"beams.{field}.{stoke.lower()}_beam{beam_num}_weight_file": filename}
    }

    return pymongo.UpdateOne(myquery, newvalues, upsert=True)


def cutout_image(
    image_name: str,
    data_in_mem: np.ndarray,
    old_header: fits.Header,
    cube: SpectralCube,
    source_id: str,
    cutout_args: CutoutArgs,
    field: str,
    beam_num: int,
    stoke: str,
    pad=3,
    dryrun=False,
) -> pymongo.UpdateOne:
    """Perform a cutout.

    Returns:
        pymongo.UpdateOne: Update query for MongoDB
    """
    logger.setLevel(logging.INFO)

    outdir = os.path.abspath(cutout_args.outdir)

    ret = []
    basename = os.path.basename(image_name)
    outname = f"{source_id}.cutout.{basename}"
    outfile = os.path.join(outdir, outname)

    padder = cube.header["BMAJ"] * u.deg * pad

    xlo = Longitude(cutout_args.ra_low * u.deg) - Longitude(padder)
    xhi = Longitude(cutout_args.ra_high * u.deg) + Longitude(padder)
    ylo = Latitude(cutout_args.dec_low * u.deg) - Latitude(padder)
    yhi = Latitude(cutout_args.dec_high * u.deg) + Latitude(padder)

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
    sub_data = data_in_mem[
        :,
        :,
        yp_lo_idx:yp_hi_idx,
        xp_lo_idx:xp_hi_idx,  # freq, Stokes, y, x
    ]
    fixed_header = fix_header(new_header, old_header)
    # Add source name to header for CASDA
    fixed_header["OBJECT"] = source_id
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
    myquery = {"Source_ID": source_id}

    filename = os.path.join(
        os.path.basename(os.path.dirname(outfile)), os.path.basename(outfile)
    )
    newvalues = {
        "$set": {f"beams.{field}.{stoke.lower()}_beam{beam_num}_image_file": filename}
    }

    return pymongo.UpdateOne(myquery, newvalues, upsert=True)


def get_args(
    comps: List[Dict],
    beam: Dict,
    island_id: str,
    outdir: str,
) -> Union[List[CutoutArgs], None]:
    """Get arguments for cutout function

    Args:
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
        List[CutoutArgs]: List of cutout arguments for cutout function
    """

    logger.setLevel(logging.INFO)

    assert beam["Source_ID"] == island_id

    if len(comps) == 0:
        logger.warning(f"Skipping island {island_id} -- no components found")
        return None

    outdir = f"{outdir}/{island_id}"
    try_mkdir(outdir)

    # Find image size
    ras: List[float] = []
    decs: List[float] = []
    majs: List[float] = []
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
        ra_high = ra_max + ra_off

        ra_min = np.min(coords.ra)
        ra_i_min = np.argmin(coords.ra)
        ra_off = Longitude(majs[ra_i_min])
        ra_low = ra_min - ra_off

        dec_max = np.max(coords.dec)
        dec_i_max = np.argmax(coords.dec)
        dec_off = Longitude(majs[dec_i_max])
        dec_high = dec_max + dec_off

        dec_min = np.min(coords.dec)
        dec_i_min = np.argmin(coords.dec)
        dec_off = Longitude(majs[dec_i_min])
        dec_low = dec_min - dec_off
    except Exception as e:
        logger.debug(f"coords are {coords=}")
        logger.debug(f"comps are {comps=}")
        raise e

    return CutoutArgs(
        ra_high=ra_high.deg,
        ra_low=ra_low.deg,
        dec_high=dec_high.deg,
        dec_low=dec_low.deg,
        outdir=outdir,
    )


def worker(
    host: str,
    epoch: int,
    beam: Dict,
    comps: List[Dict],
    outdir: str,
    image_name: str,
    data_in_mem: np.ndarray,
    old_header: fits.Header,
    cube: SpectralCube,
    field: str,
    beam_num: int,
    stoke: str,
    pad: float = 3,
    username: Optional[str] = None,
    password: Optional[str] = None,
):
    _, _, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )
    cut_args = get_args(
        comps=comps,
        beam=beam,
        island_id=beam["Source_ID"],
        outdir=outdir,
    )
    image_update = cutout_image(
        image_name=image_name,
        data_in_mem=data_in_mem,
        old_header=old_header,
        cube=cube,
        source_id=beam["Source_ID"],
        cutout_args=cut_args,
        field=field,
        beam_num=beam_num,
        stoke=stoke,
        pad=pad,
        dryrun=False,
    )
    weight_update = cutout_weight(
        image_name=image_name,
        source_id=beam["Source_ID"],
        cutout_args=cut_args,
        field=field,
        beam_num=beam_num,
        stoke=stoke,
        dryrun=False,
    )
    return [image_update, weight_update]


@task(name="Cutout from big cube")
def big_cutout(
    beams: List[Dict],
    beam_num: int,
    stoke: str,
    datadir: str,
    outdir: str,
    host: str,
    epoch: int,
    field: str,
    pad: float = 3,
    username: Optional[str] = None,
    password: Optional[str] = None,
    limit: Optional[int] = None,
) -> List[pymongo.UpdateOne]:
    with open("comps.pkl", "rb") as f:
        comps_dict = pickle.load(f)
    wild = (
        f"{datadir}/image.restored.{stoke.lower()}*contcube*beam{beam_num:02}.conv.fits"
    )
    images = glob(wild)
    if len(images) == 0:
        raise Exception(f"No images found matching '{wild}'")
    elif len(images) > 1:
        raise Exception(f"More than one image found matching '{wild}'. Files {images=}")

    image_name = images[0]

    # Read the whole lad into memory
    logger.info(f"Reading {image_name}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        cube = SpectralCube.read(image_name, memmap=True, mode="denywrite")

    data_in_mem = np.array(fits.getdata(image_name))
    old_header = fits.getheader(image_name)

    if limit is not None:
        logger.critical(f"Limiting to {limit} islands")
        beams = beams[:limit]

    updates: List[pymongo.UpdateOne] = []
    with ThreadPoolExecutor() as executor:
        futures = []
        for beam in beams:
            futures.append(
                executor.submit(
                    worker,
                    host=host,
                    epoch=epoch,
                    beam=beam,
                    comps=comps_dict[beam["Source_ID"]],
                    outdir=outdir,
                    image_name=image_name,
                    data_in_mem=data_in_mem,
                    old_header=old_header,
                    cube=cube,
                    field=field,
                    beam_num=beam_num,
                    stoke=stoke,
                    pad=pad,
                    username=username,
                    password=password,
                )
            )
        for future in tqdm(futures, file=TQDM_OUT, desc=f"Cutting {image_name}"):
            updates += future.result()

    return updates


@flow(name="Cutout islands")
def cutout_islands(
    field: str,
    directory: str,
    host: str,
    epoch: int,
    sbid: Optional[int] = None,
    username: Optional[str] = None,
    password: Optional[str] = None,
    pad: float = 3,
    stokeslist: Optional[List[str]] = None,
    dryrun: bool = True,
    limit: Optional[int] = None,
) -> None:
    """Flow to cutout islands in parallel.

    Args:
        field (str): RACS field name.
        directory (str): Directory to store cutouts.
        host (str): MongoDB host.
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo password. Defaults to None.
        verbose (bool, optional): Verbose output. Defaults to True.
        pad (int, optional): Number of beamwidths to pad cutouts. Defaults to 3.
        stokeslist (List[str], optional): Stokes parameters to cutout. Defaults to None.
        dryrun (bool, optional): Do everything except write FITS files. Defaults to True.
    """
    if stokeslist is None:
        stokeslist = ["I", "Q", "U", "V"]

    directory = os.path.abspath(directory)
    outdir = os.path.join(directory, "cutouts")

    logger.info("Testing database. ")
    test_db(
        host=host,
        username=username,
        password=password,
    )

    beams_col, island_col, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )

    # Query the DB

    query = {"$and": [{f"beams.{field}": {"$exists": True}}]}
    if sbid is not None:
        query["$and"].append({f"beams.{field}.SBIDs": sbid})

    unique_beams_nums: Set[int] = set(
        beams_col.distinct(f"beams.{field}.beam_list", query)
    )
    source_ids = sorted(beams_col.distinct("Source_ID", query))

    beams_dict: Dict[int, List[Dict]] = {b: [] for b in unique_beams_nums}

    query = {
        "$and": [
            {f"beams.{field}": {"$exists": True}},
            {f"beams.{field}.beam_list": {"$in": list(unique_beams_nums)}},
        ]
    }
    if sbid is not None:
        query["$and"].append({f"beams.{field}.SBIDs": sbid})

    all_beams = list(beams_col.find(query).sort("Source_ID"))
    for beams in tqdm(all_beams, desc="Getting beams", file=TQDM_OUT):
        for beam_num in beams[f"beams"][field]["beam_list"]:
            beams_dict[beam_num].append(beams)

    comps_dict: Dict[str, List[Dict]] = {s: [] for s in source_ids}
    all_comps = list(
        comp_col.find({"Source_ID": {"$in": source_ids}}).sort("Source_ID")
    )
    for comp in tqdm(all_comps, desc="Getting components", file=TQDM_OUT):
        comps_dict[comp["Source_ID"]].append(comp)

    # Dump comps to file
    with open("comps.pkl", "wb") as f:
        pickle.dump(comps_dict, f)

    # Create output dir if it doesn't exist
    try_mkdir(outdir)
    cuts: List[pymongo.UpdateOne] = []
    for stoke in stokeslist:
        for beam_num in unique_beams_nums:
            results = big_cutout.submit(
                beams=beams_dict[beam_num],
                beam_num=beam_num,
                stoke=stoke,
                datadir=directory,
                outdir=outdir,
                host=host,
                epoch=epoch,
                field=field,
                pad=pad,
                username=username,
                password=password,
                limit=limit,
            )
            cuts.append(results)

    if not dryrun:
        _updates = [f.result() for f in cuts]
        updates = [val for sublist in _updates for val in sublist]
        logger.info("Updating database...")
        db_res = beams_col.bulk_write(updates, ordered=False)
        logger.info(pformat(db_res.bulk_api_result))

    os.remove("comps.pkl")

    logger.info("Cutouts Done!")


def main(args: argparse.Namespace) -> None:
    """Main script

    Args:
        args (argparse.Namespace): Command-line args
        verbose (bool, optional): Verbose output. Defaults to True.
    """
    cutout_islands(
        field=args.field,
        directory=args.datadir,
        host=args.host,
        epoch=args.epoch,
        sbid=args.sbid,
        username=args.username,
        password=args.password,
        pad=args.pad,
        stokeslist=args.stokeslist,
        dryrun=args.dryrun,
        limit=args.limit,
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
        formatter_class=UltimateHelpFormatter,
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
        "-e",
        "--epoch",
        type=int,
        default=0,
        help="Epoch of observation.",
    )

    parser.add_argument(
        "--sbid",
        type=int,
        default=None,
        help="SBID of observation.",
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
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit number of islands to process [None]",
    )

    return cut_parser


def cli() -> None:
    """Command-line interface"""
    parser = cutout_parser()

    args = parser.parse_args()

    verbose = args.verbose
    if verbose:
        logger.setLevel(logging.INFO)

    test_db(
        host=args.host,
        username=args.username,
        password=args.password,
    )

    main(args)


if __name__ == "__main__":
    cli()
