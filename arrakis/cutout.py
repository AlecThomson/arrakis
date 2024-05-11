#!/usr/bin/env python
"""Produce cutouts from RACS cubes"""

import argparse
import logging
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from pprint import pformat
from shutil import copyfile
from typing import List
from typing import NamedTuple as Struct
from typing import Optional, Set, TypeVar

import astropy.units as u
import numpy as np
import pandas as pd
import pymongo
from astropy.coordinates import Latitude, Longitude, SkyCoord
from astropy.io import fits
from astropy.utils import iers
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs.utils import skycoord_to_pixel
from prefect import flow, task
from spectral_cube import SpectralCube
from spectral_cube.utils import SpectralCubeWarning
from tqdm.auto import tqdm

from arrakis.logger import TqdmToLogger, UltimateHelpFormatter, logger
from arrakis.utils.database import (
    get_db,
    get_field_db,
    test_db,
    validate_sbid_field_pair,
)
from arrakis.utils.fitsutils import fix_header
from arrakis.utils.pipeline import generic_parser, logo_str, workdir_arg_parser

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
    outdir: Path
    """Output directory"""


def cutout_weight(
    image_name: Path,
    source_id: str,
    cutout_args: Optional[CutoutArgs],
    field: str,
    stoke: str,
    beam_num: int,
    dryrun=False,
) -> pymongo.UpdateOne:
    # Update database
    myquery = {"Source_ID": source_id}

    if cutout_args is None:
        logger.error(f"Skipping {source_id} -- no components found")
        newvalues = {
            "$set": {f"beams.{field}.{stoke.lower()}_beam{beam_num}_weight_file": ""}
        }
        return pymongo.UpdateOne(myquery, newvalues, upsert=True)

    outdir = cutout_args.outdir.absolute()
    basename = image_name.name
    outname = f"{source_id}.cutout.{basename}"
    outfile = outdir / outname
    image = (
        image_name.parent
        / image_name.name.replace("image.restored", "weights.restored")
    ).with_suffix(".txt")
    outfile = (
        outfile.parent / outfile.name.replace("image.restored", "weights.restored")
    ).with_suffix(".txt")

    if not dryrun:
        copyfile(image, outfile)
        logger.info(f"Written to {outfile}")

    filename = outfile.parent / outfile.name
    newvalues = {
        "$set": {
            f"beams.{field}.{stoke.lower()}_beam{beam_num}_weight_file": filename.absolute().as_posix()
        }
    }

    return pymongo.UpdateOne(myquery, newvalues, upsert=True)


def cutout_image(
    image_name: Path,
    data_in_mem: np.ndarray,
    old_header: fits.Header,
    cube: SpectralCube,
    source_id: str,
    cutout_args: Optional[CutoutArgs],
    field: str,
    beam_num: int,
    stoke: str,
    pad: float = 3,
    dryrun: bool = False,
) -> pymongo.UpdateOne:
    """Perform a cutout.

    Returns:
        pymongo.UpdateOne: Update query for MongoDB
    """
    logger.setLevel(logging.INFO)
    # Update database
    myquery = {"Source_ID": source_id}
    if cutout_args is None:
        logger.error(f"Skipping {source_id} -- no components found")
        newvalues = {
            "$set": {f"beams.{field}.{stoke.lower()}_beam{beam_num}_weight_file": ""}
        }
        return pymongo.UpdateOne(myquery, newvalues, upsert=True)

    outdir = cutout_args.outdir.absolute()

    basename = image_name.name
    outname = f"{source_id}.cutout.{basename}"
    outfile = outdir / outname

    padder: float = cube.header["BMAJ"] * u.deg * pad

    top_right = SkyCoord(cutout_args.ra_low * u.deg, cutout_args.dec_high * u.deg)
    bottom_left = SkyCoord(cutout_args.ra_high * u.deg, cutout_args.dec_low * u.deg)

    top_right_off = top_right.spherical_offsets_by(d_lon=padder, d_lat=padder)
    bottom_left_off = bottom_left.spherical_offsets_by(d_lon=-padder, d_lat=-padder)
    # Only need the critical corners - but just in case:
    # top_left = SkyCoord(cutout_args.ra_high * u.deg, cutout_args.dec_high * u.deg)
    # bottom_right = SkyCoord(cutout_args.ra_low, cutout_args.dec_low * u.deg)
    # top_left_off = top_left.spherical_offsets_by(d_lon=-padder, d_lat=padder)
    # bottom_right_off = bottom_right.spherical_offsets_by(d_lon=padder, d_lat=-padder)

    x_left, y_bottom = skycoord_to_pixel(bottom_left_off, cube.wcs)
    x_right, y_top = skycoord_to_pixel(top_right_off, cube.wcs)

    # Round for cutout
    yp_lo_idx = int(np.floor(y_bottom))
    yp_hi_idx = int(np.ceil(y_top))
    xp_lo_idx = int(np.floor(x_right))
    xp_hi_idx = int(np.ceil(x_left))

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

    filename = outfile.parent / outfile.name
    newvalues = {
        "$set": {
            f"beams.{field}.{stoke.lower()}_beam{beam_num}_image_file": filename.absolute().as_posix()
        }
    }

    return pymongo.UpdateOne(myquery, newvalues, upsert=True)


def get_args(
    comps: pd.DataFrame,
    source: pd.Series,
    outdir: Path,
) -> Optional[CutoutArgs]:
    """Get arguments for cutout function

    Args:
        comps (pd.DataFrame): List of mongo entries for RACS components in island
        beam (Dict): Mongo entry for the RACS beam
        island_id (str): RACS island ID
        outdir (Path): Input directory

    Raises:
        e: Exception
        Exception: Problems with coordinates

    Returns:
        List[CutoutArgs]: List of cutout arguments for cutout function
    """

    logger.setLevel(logging.INFO)

    island_id = source.Source_ID

    if len(comps) == 0:
        logger.warning(f"Skipping island {island_id} -- no components found")
        return None

    outdir = outdir / island_id
    outdir.mkdir(parents=True, exist_ok=True)

    # Find image size
    ras: u.Quantity = comps.RA.values * u.deg
    decs: u.Quantity = comps.Dec.values * u.deg
    majs: List[float] = comps.Maj.values * u.arcsec

    coords = SkyCoord(ras, decs)

    try:
        ra_max = np.max(coords.ra)
        ra_i_max = np.argmax(coords.ra)
        ra_off = Longitude(majs[ra_i_max])

        ra_min = np.min(coords.ra)
        ra_i_min = np.argmin(coords.ra)
        ra_off = Longitude(majs[ra_i_min])

        dec_max = np.max(coords.dec)
        dec_i_max = np.argmax(coords.dec)
        dec_off = Latitude(majs[dec_i_max])

        dec_min = np.min(coords.dec)
        dec_i_min = np.argmin(coords.dec)
        dec_off = Latitude(majs[dec_i_min])

        # Use SkyCoords to account for poles and meridian
        top_right = SkyCoord(ra_min, dec_max)
        top_left = SkyCoord(ra_max, dec_max)
        bottom_right = SkyCoord(ra_min, dec_min)
        bottom_left = SkyCoord(ra_max, dec_min)
        top_right_off = top_right.spherical_offsets_by(d_lon=ra_off, d_lat=dec_off)
        top_left_off = top_left.spherical_offsets_by(d_lon=-ra_off, d_lat=dec_off)
        bottom_right_off = bottom_right.spherical_offsets_by(
            d_lon=ra_off, d_lat=-dec_off
        )
        bottom_left_off = bottom_left.spherical_offsets_by(
            d_lon=-ra_off, d_lat=-dec_off
        )
        ra_high: float = np.max(
            [
                top_right_off.ra.deg,
                top_left_off.ra.deg,
                bottom_right_off.ra.deg,
                bottom_left_off.ra.deg,
            ]
        )
        ra_low: float = np.min(
            [
                top_right_off.ra.deg,
                top_left_off.ra.deg,
                bottom_right_off.ra.deg,
                bottom_left_off.ra.deg,
            ]
        )
        dec_high: float = np.max(
            [
                top_right_off.dec.deg,
                top_left_off.dec.deg,
                bottom_right_off.dec.deg,
                bottom_left_off.dec.deg,
            ]
        )
        dec_low: float = np.min(
            [
                top_right_off.dec.deg,
                top_left_off.dec.deg,
                bottom_right_off.dec.deg,
                bottom_left_off.dec.deg,
            ]
        )

    except Exception as e:
        logger.debug(f"coords are {coords=}")
        logger.debug(f"comps are {comps=}")
        raise e

    return CutoutArgs(
        ra_high=ra_high,
        ra_low=ra_low,
        dec_high=dec_high,
        dec_low=dec_low,
        outdir=outdir,
    )


def worker(
    host: str,
    epoch: int,
    source: pd.Series,
    comps: pd.DataFrame,
    outdir: Path,
    image_name: Path,
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
        source=source,
        outdir=outdir,
    )
    image_update = cutout_image(
        image_name=image_name,
        data_in_mem=data_in_mem,
        old_header=old_header,
        cube=cube,
        source_id=source.Source_ID,
        cutout_args=cut_args,
        field=field,
        beam_num=beam_num,
        stoke=stoke,
        pad=pad,
        dryrun=False,
    )
    weight_update = cutout_weight(
        image_name=image_name,
        source_id=source.Source_ID,
        cutout_args=cut_args,
        field=field,
        beam_num=beam_num,
        stoke=stoke,
        dryrun=False,
    )
    return [image_update, weight_update]


@task(name="Cutout from big cube")
def big_cutout(
    sources: pd.DataFrame,
    comps: pd.DataFrame,
    beam_num: int,
    stoke: str,
    datadir: Path,
    outdir: Path,
    host: str,
    epoch: int,
    field: str,
    pad: float = 3,
    username: Optional[str] = None,
    password: Optional[str] = None,
    limit: Optional[int] = None,
) -> List[pymongo.UpdateOne]:
    wild = f"image.restored.{stoke.lower()}*contcube*beam{beam_num:02}.conv.fits"
    images = list(datadir.glob(wild))
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
        sources = sources[:limit]

    updates: List[pymongo.UpdateOne] = []
    with ThreadPoolExecutor() as executor:
        futures = []
        for _, source in sources.iterrows():
            futures.append(
                executor.submit(
                    worker,
                    host=host,
                    epoch=epoch,
                    source=source,
                    comps=comps.loc[source],
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
    directory: Path,
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
        directory (Path): Directory to store cutouts.
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

    directory = directory.absolute()
    outdir = directory / "cutouts"

    logger.info("Testing database. ")
    test_db(
        host=host,
        username=username,
        password=password,
    )

    beams_col, island_col, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )

    # Check for SBID match
    if sbid is not None:
        field_col = get_field_db(
            host=host,
            epoch=epoch,
            username=username,
            password=password,
        )
        sbid_check = validate_sbid_field_pair(
            field_name=field,
            sbid=sbid,
            field_col=field_col,
        )
        if not sbid_check:
            raise ValueError(f"SBID {sbid} does not match field {field}")

    query = {"$and": [{f"beams.{field}": {"$exists": True}}]}
    if sbid is not None:
        query["$and"].append({f"beams.{field}.SBIDs": sbid})

    unique_beams_nums: Set[int] = set(
        beams_col.distinct(f"beams.{field}.beam_list", query)
    )
    source_ids = sorted(beams_col.distinct("Source_ID", query))

    # beams_dict: Dict[int, List[Dict]] = {b: [] for b in unique_beams_nums}

    query = {
        "$and": [
            {f"beams.{field}": {"$exists": True}},
            {f"beams.{field}.beam_list": {"$in": list(unique_beams_nums)}},
        ]
    }
    if sbid is not None:
        query["$and"].append({f"beams.{field}.SBIDs": sbid})

    beams_df = pd.DataFrame(
        beams_col.find(query, {"Source_ID": 1, f"beams.{field}.beam_list": 1}).sort(
            "Source_ID"
        )
    )

    beam_source_list = []
    for i, row in tqdm(beams_df.iterrows()):
        beam_list = row.beams[field]["beam_list"]
        for b in beam_list:
            beam_source_list.append({"Source_ID": row.Source_ID, "beam": b})
    beam_source_df = pd.DataFrame(beam_source_list)
    beam_source_df.set_index("beam", inplace=True)

    comps_df = pd.DataFrame(
        comp_col.find({"Source_ID": {"$in": source_ids}}).sort("Source_ID")
    )
    comps_df.set_index("Source_ID", inplace=True)

    # Create output dir if it doesn't exist
    outdir.mkdir(parents=True, exist_ok=True)
    cuts: List[pymongo.UpdateOne] = []
    for stoke in stokeslist:
        for beam_num in unique_beams_nums:
            results = big_cutout.submit(
                sources=beam_source_df.loc[beam_num],
                comps=comps_df.loc[beam_source_df.loc[beam_num].Source_ID],
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

    return cut_parser


def cli() -> None:
    """Command-line interface"""
    gen_parser = generic_parser(parent_parser=True)
    work_parser = workdir_arg_parser(parent_parser=True)
    cut_parser = cutout_parser(parent_parser=True)
    parser = argparse.ArgumentParser(
        formatter_class=UltimateHelpFormatter,
        parents=[gen_parser, work_parser, cut_parser],
        description=cut_parser.description,
    )
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
