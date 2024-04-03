#!/usr/bin/env python3
"""Merge multiple RACS fields"""
import os
from pprint import pformat
from shutil import copyfile
from typing import Dict, List, Optional

import pymongo
from prefect import flow, task, unmapped

from arrakis.linmos import get_yanda, linmos
from arrakis.logger import UltimateHelpFormatter, logger
from arrakis.utils.database import get_db, test_db
from arrakis.utils.io import try_mkdir


def make_short_name(name: str) -> str:
    short = os.path.join(
        os.path.basename(os.path.dirname(name)), os.path.basename(name)
    )
    return short


@task(name="Copy singleton island")
def copy_singleton(
    beam: dict, field_dict: Dict[str, str], merge_name: str, data_dir: str
) -> List[pymongo.UpdateOne]:
    """Copy an island within a single field to the merged field

    Args:
        beam (dict): Beam document
        field_dict (Dict[str, str]): Field dictionary
        merge_name (str): Merged field name
        data_dir (str): Output directory

    Raises:
        KeyError: If ion files not found

    Returns:
        List[pymongo.UpdateOne]: Database updates
    """
    updates = []
    for field, vals in beam["beams"].items():
        if field not in field_dict.keys():
            continue
        field_dir = field_dict[field]
        try:
            i_file_old = os.path.join(field_dir, vals["i_file"])
            q_file_old = os.path.join(field_dir, vals["q_file_ion"])
            u_file_old = os.path.join(field_dir, vals["u_file_ion"])
        except KeyError:
            raise KeyError("Ion files not found. Have you run FRion?")
        new_dir = os.path.join(data_dir, beam["Source_ID"])

        try_mkdir(new_dir, verbose=False)

        i_file_new = os.path.join(new_dir, os.path.basename(i_file_old)).replace(
            ".fits", ".edge.linmos.fits"
        )
        q_file_new = os.path.join(new_dir, os.path.basename(q_file_old)).replace(
            ".fits", ".edge.linmos.fits"
        )
        u_file_new = os.path.join(new_dir, os.path.basename(u_file_old)).replace(
            ".fits", ".edge.linmos.fits"
        )

        for src, dst in zip(
            [i_file_old, q_file_old, u_file_old], [i_file_new, q_file_new, u_file_new]
        ):
            copyfile(src, dst)
            src_weight = src.replace(".image.restored.", ".weights.").replace(
                ".ion", ""
            )
            dst_weight = dst.replace(".image.restored.", ".weights.").replace(
                ".ion", ""
            )
            copyfile(src_weight, dst_weight)

        query = {"Source_ID": beam["Source_ID"]}
        newvalues = {
            "$set": {
                f"beams.{merge_name}.i_file": make_short_name(i_file_new),
                f"beams.{merge_name}.q_file": make_short_name(q_file_new),
                f"beams.{merge_name}.u_file": make_short_name(u_file_new),
                f"beams.{merge_name}.DR1": True,
            }
        }

        updates.append(pymongo.UpdateOne(query, newvalues))
    return updates


def copy_singletons(
    field_dict: Dict[str, str],
    data_dir: str,
    beams_col: pymongo.collection.Collection,
    merge_name: str,
) -> List[pymongo.UpdateOne]:
    """Copy islands that don't overlap other fields

    Args:
        field_dict (Dict[str, str]): Field dictionary
        data_dir (str): Data directory
        beams_col (pymongo.collection.Collection): Beams collection
        merge_name (str): Merged field name

    Returns:
        List[pymongo.UpdateOne]: Database updates
    """
    # Find all islands with the given fields that DON'T overlap another field
    query = {
        "$or": [
            {
                "$and": [
                    {f"beams.{field}": {"$exists": True}},
                    {f"beams.{field}.DR1": True},
                    {"n_fields_DR1": 1},
                ]
            }
            for field in field_dict.keys()
        ]
    }

    island_ids = sorted(beams_col.distinct("Source_ID", query))
    big_beams = list(
        beams_col.find({"Source_ID": {"$in": island_ids}}).sort("Source_ID")
    )
    updates = copy_singleton.map(
        beam=big_beams,
        field_dict=unmapped(field_dict),
        merge_name=unmapped(merge_name),
        data_dir=unmapped(data_dir),
    )
    return updates


def genparset(
    old_ims: list,
    stokes: str,
    new_dir: str,
) -> str:
    imlist = "[" + ",".join([im.replace(".fits", "") for im in old_ims]) + "]"
    weightlist = f"[{','.join([im.replace('.fits', '').replace('.image.restored.','.weights.').replace('.ion','') for im in old_ims])}]"

    im_outname = os.path.join(new_dir, os.path.basename(old_ims[0])).replace(
        ".fits", ".edge.linmos"
    )
    wt_outname = (
        os.path.join(new_dir, os.path.basename(old_ims[0]))
        .replace(".fits", ".edge.linmos")
        .replace(".image.restored.", ".weights.")
    )

    parset_file = os.path.join(new_dir, f"edge_linmos_{stokes}.in")
    parset = f"""# LINMOS parset
linmos.names            = {imlist}
linmos.weights          = {weightlist}
linmos.imagetype        = fits
linmos.outname          = {im_outname}
linmos.outweight        = {wt_outname}
# For ASKAPsoft>1.3.0
linmos.weighttype       = FromWeightImages
linmos.weightstate      = Corrected
"""

    with open(parset_file, "w") as f:
        f.write(parset)

    return parset_file


def merge_multiple_field(
    beam: dict, field_dict: dict, merge_name: str, data_dir: str, image: str
) -> List[pymongo.UpdateOne]:
    """Merge an island that overlaps multiple fields

    Args:
        beam (dict): Beam document
        field_dict (dict): Field dictionary
        merge_name (str): Merged field name
        data_dir (str): Data directory
        image (str): Yandasoft image

    Raises:
        KeyError: If ion files not found

    Returns:
        List[pymongo.UpdateOne]: Database updates
    """
    i_files_old = []
    q_files_old = []
    u_files_old = []
    updates = []
    for field, vals in beam["beams"].items():
        if field not in field_dict.keys():
            continue
        field_dir = field_dict[field]
        try:
            i_file_old = os.path.join(field_dir, vals["i_file"])
            q_file_old = os.path.join(field_dir, vals["q_file_ion"])
            u_file_old = os.path.join(field_dir, vals["u_file_ion"])
        except KeyError:
            raise KeyError("Ion files not found. Have you run FRion?")
        i_files_old.append(i_file_old)
        q_files_old.append(q_file_old)
        u_files_old.append(u_file_old)

    new_dir = os.path.join(data_dir, beam["Source_ID"])

    try_mkdir(new_dir, verbose=False)

    for stokes, imlist in zip(["I", "Q", "U"], [i_files_old, q_files_old, u_files_old]):
        parset_file = genparset(imlist, stokes, new_dir)
        update = linmos.fn(parset_file, merge_name, image)
        updates.append(update)

    return updates


@task(name="Merge multiple islands")
def merge_multiple_fields(
    field_dict: Dict[str, str],
    data_dir: str,
    beams_col: pymongo.collection.Collection,
    merge_name: str,
    image: str,
) -> List[pymongo.UpdateOne]:
    """Merge multiple islands that overlap multiple fields

    Args:
        field_dict (Dict[str, str]): Field dictionary
        data_dir (str): Data directory
        beams_col (pymongo.collection.Collection): Beams collection
        merge_name (str): Merged field name
        image (str): Yandasoft image

    Returns:
        List[pymongo.UpdateOne]: Database updates
    """
    # Find all islands with the given fields that overlap another field
    query = {
        "$or": [
            {
                "$and": [
                    {f"beams.{field}": {"$exists": True}},
                    {f"beams.{field}.DR1": True},
                    {"n_fields_DR1": {"$gt": 1}},
                ]
            }
            for field in field_dict.keys()
        ]
    }

    island_ids = sorted(beams_col.distinct("Source_ID", query))
    big_beams = list(
        beams_col.find({"Source_ID": {"$in": island_ids}}).sort("Source_ID")
    )

    updates = merge_multiple_field.map(
        beam=big_beams,
        field_dict=unmapped(field_dict),
        merge_name=unmapped(merge_name),
        data_dir=unmapped(data_dir),
        image=unmapped(image),
    )

    return updates


@flow(name="Merge fields")
def main(
    fields: List[str],
    field_dirs: List[str],
    merge_name: str,
    output_dir: str,
    host: str,
    epoch: int,
    username: Optional[str] = None,
    password: Optional[str] = None,
    yanda="1.3.0",
) -> str:
    logger.debug(f"{fields=}")

    assert len(fields) == len(
        field_dirs
    ), f"List of fields must be the same length as length of field dirs. {len(fields)=},{len(field_dirs)=}"

    field_dict = {
        field: os.path.join(field_dir, "cutouts")
        for field, field_dir in zip(fields, field_dirs)
    }

    image = get_yanda(version=yanda)

    beams_col, island_col, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )

    output_dir = os.path.abspath(output_dir)
    inter_dir = os.path.join(output_dir, merge_name)
    try_mkdir(inter_dir)
    data_dir = os.path.join(inter_dir, "cutouts")
    try_mkdir(data_dir)

    singleton_updates = copy_singletons(
        field_dict=field_dict,
        data_dir=data_dir,
        beams_col=beams_col,
        merge_name=merge_name,
    )

    mutilple_updates = merge_multiple_fields(
        field_dict=field_dict,
        data_dir=data_dir,
        beams_col=beams_col,
        merge_name=merge_name,
        image=image,
    )

    singleton_comp = [f.result() for f in singleton_updates]
    multiple_comp = [f.result() for f in mutilple_updates]

    for m in multiple_comp:
        m._doc["$set"].update({f"beams.{merge_name}.DR1": True})

    db_res_single = beams_col.bulk_write(singleton_comp, ordered=False)
    logger.info(pformat(db_res_single.bulk_api_result))

    db_res_multiple = beams_col.bulk_write(multiple_comp, ordered=False)
    logger.info(pformat(db_res_multiple.bulk_api_result))

    logger.info("LINMOS Done!")
    return inter_dir


def cli():
    """Command-line interface"""
    import argparse

    # Help string to be shown using the -h option
    descStr = """
    Mosaic RACS beam fields with linmos.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=UltimateHelpFormatter
    )

    parser.add_argument(
        "--merge_name",
        type=str,
        help="Name of the merged region",
    )

    parser.add_argument(
        "--fields", type=str, nargs="+", help="RACS fields to mosaic - e.g. 2132-50A."
    )

    parser.add_argument(
        "--datadirs",
        type=str,
        nargs="+",
        help="Directories containing cutouts (in subdir outdir/cutouts)..",
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path to save merged data (in output_dir/merge_name/cutouts)",
    )

    parser.add_argument(
        "--yanda",
        type=str,
        default="1.3.0",
        help="Yandasoft version to pull from DockerHub [1.3.0].",
    )

    parser.add_argument(
        "--host",
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
        "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )

    args = parser.parse_args()
    verbose = args.verbose
    test_db(
        host=args.host, username=args.username, password=args.password, verbose=verbose
    )

    main(
        fields=args.fields,
        field_dirs=args.datadirs,
        merge_name=args.merge_name,
        output_dir=args.output_dir,
        host=args.host,
        epoch=args.epoch,
        username=args.username,
        password=args.password,
        yanda=args.yanda,
    )


if __name__ == "__main__":
    cli()
