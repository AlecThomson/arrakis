#!/usr/bin/env python3
"""Merge multiple RACS fields."""

from __future__ import annotations

import argparse
from pathlib import Path
from pprint import pformat
from shutil import copyfile
from typing import Any

import pymongo
from prefect import flow, task

from arrakis.linmos import get_yanda, linmos, linmos_parser
from arrakis.logger import UltimateHelpFormatter, logger
from arrakis.utils.database import get_db, test_db
from arrakis.utils.io import try_mkdir


def make_short_name(name: Path) -> str:
    """Make a short name for a file.

    Args:
        name (Path): File name

    Returns:
        str: Short name
    """
    return (Path(name.parent.name) / name.name).as_posix()


@task(name="Copy singleton island")
def copy_singleton(
    beam: dict[str, Any], field_dict: dict[str, Path], merge_name: str, data_dir: Path
) -> list[pymongo.UpdateOne]:
    """Copy an island within a single field to the merged field.

    Args:
        beam (dict): Beam document
        field_dict (dict[str, Path]): Field dictionary
        merge_name (str): Merged field name
        data_dir (Path): Output directory

    Raises:
        KeyError: If ion files not found

    Returns:
        List[pymongo.UpdateOne]: Database updates
    """
    updates = []
    for field, vals in beam["beams"].items():
        if field not in field_dict:
            continue
        field_dir = field_dict[field]
        try:
            i_file_old = field_dir / str(vals["i_file"])
            q_file_old = field_dir / str(vals["q_file_ion"])
            u_file_old = field_dir / str(vals["u_file_ion"])
        except KeyError as e:
            msg = "Ion files not found. Have you run FRion?"
            raise KeyError(msg) from e
        new_dir = data_dir / str(beam["Source_ID"])
        new_dir.mkdir(exist_ok=True)

        i_file_new = (new_dir / i_file_old.name).with_suffix(".edge.linmos.fits")
        q_file_new = (new_dir / q_file_old.name).with_suffix(".edge.linmos.fits")
        u_file_new = (new_dir / u_file_old.name).with_suffix(".edge.linmos.fits")

        for src, dst in zip(
            [i_file_old, q_file_old, u_file_old], [i_file_new, q_file_new, u_file_new]
        ):
            copyfile(src, dst)
            src_weight = (
                src.as_posix()
                .replace(".image.restored.", ".weights.")
                .replace(".ion", "")
            )
            dst_weight = (
                dst.as_posix()
                .replace(".image.restored.", ".weights.")
                .replace(".ion", "")
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
    field_dict: dict[str, str],
    data_dir: str,
    beams_col: pymongo.collection.Collection,
    merge_name: str,
) -> list[pymongo.UpdateOne]:
    """Copy islands that don't overlap other fields.

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
            for field in field_dict
        ]
    }

    island_ids = sorted(beams_col.distinct("Source_ID", query))
    big_beams = list(
        beams_col.find({"Source_ID": {"$in": island_ids}}).sort("Source_ID")
    )
    updates = []
    for beam in big_beams:
        update = copy_singleton.submit(
            beam=beam,
            field_dict=field_dict,
            merge_name=merge_name,
            data_dir=data_dir,
        )
        updates.append(update)
    return updates


def genparset(
    old_ims: list[Path],
    stokes: str,
    new_dir: Path,
) -> Path:
    """Generate a linmos parset file.

    Args:
        old_ims (list[Path]): Old images
        stokes (str): Stokes parameter
        new_dir (Path): Output directory

    Returns:
        Path: Path to parset file
    """
    imlist = "[" + ",".join([im.with_suffix("").as_posix() for im in old_ims]) + "]"
    weightlist = f"[{','.join([im.with_suffix("").as_posix().replace('.image.restored.','.weights.').replace('.ion','') for im in old_ims])}]"

    im_outname = (new_dir / old_ims[0].name).with_suffix(".edge.linmos").as_posix()
    wt_outname = (
        (new_dir / old_ims[0].name)
        .with_suffix(".edge.linmos")
        .as_posix()
        .replace(".image.restored.", ".weights.")
    )

    parset_file = new_dir / f"edge_linmos_{stokes}.in"
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

    with parset_file.open("w") as f:
        f.write(parset)

    return parset_file


def merge_multiple_field(
    beam: dict, field_dict: dict[str, Path], merge_name: str, data_dir: Path, image: str
) -> list[pymongo.UpdateOne]:
    """Merge an island that overlaps multiple fields.

    Args:
        beam (dict): Beam document
        field_dict (dict): Field dictionary
        merge_name (str): Merged field name
        data_dir (Path): Data directory
        image (str): Yandasoft image

    Raises:
        KeyError: If ion files not found

    Returns:
        List[pymongo.UpdateOne]: Database updates
    """
    i_files_old: list[Path] = []
    q_files_old: list[Path] = []
    u_files_old: list[Path] = []
    updates = []
    for field, vals in beam["beams"].items():
        if field not in field_dict:
            continue
        field_dir = field_dict[field]
        try:
            i_file_old = field_dir / str(vals["i_file"])
            q_file_old = field_dir / str(vals["q_file_ion"])
            u_file_old = field_dir / str(vals["u_file_ion"])
        except KeyError as e:
            msg = "Ion files not found. Have you run FRion?"
            raise KeyError(msg) from e
        i_files_old.append(i_file_old)
        q_files_old.append(q_file_old)
        u_files_old.append(u_file_old)

    new_dir = data_dir / beam["Source_ID"]

    try_mkdir(new_dir, verbose=False)

    for stokes, imlist in zip(["I", "Q", "U"], [i_files_old, q_files_old, u_files_old]):
        parset_file = genparset(imlist, stokes, new_dir)
        update = linmos.fn(parset_file, merge_name, image)
        updates.append(update)

    return updates


@task(name="Merge multiple islands")
def merge_multiple_fields(
    field_dict: dict[str, str],
    data_dir: str,
    beams_col: pymongo.collection.Collection,
    merge_name: str,
    image: str,
) -> list[pymongo.UpdateOne]:
    """Merge multiple islands that overlap multiple fields.

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
            for field in field_dict
        ]
    }

    island_ids = sorted(beams_col.distinct("Source_ID", query))
    big_beams = list(
        beams_col.find({"Source_ID": {"$in": island_ids}}).sort("Source_ID")
    )
    updates = []
    for beam in big_beams:
        update = merge_multiple_field.submit(
            beam=beam,
            field_dict=field_dict,
            merge_name=merge_name,
            data_dir=data_dir,
            image=image,
        )
        updates.append(update)

    return updates


@flow(name="Merge fields")
def main(
    fields: list[str],
    field_dirs: list[Path],
    merge_name: str,
    output_dir: Path,
    host: str,
    epoch: int,
    username: str | None = None,
    password: str | None = None,
    yanda="1.3.0",
) -> str:
    """Merge multiple RACS fields.

    Args:
        fields (list[str]): List of field names.
        field_dirs (list[Path]): List of field directories.
        merge_name (str): Name of the merged field.
        output_dir (Path): Output directory.
        host (str): MongoDB host.
        epoch (int): Epoch.
        username (str | None, optional): MongoDB username. Defaults to None.
        password (str | None, optional): MongoDB password. Defaults to None.
        yanda (str, optional): Yandasoft version. Defaults to "1.3.0".

    Returns:
        str: Intermediate directory
    """
    logger.debug(f"{fields=}")

    assert (
        len(fields) == len(field_dirs)
    ), f"List of fields must be the same length as length of field dirs. {len(fields)=},{len(field_dirs)=}"

    field_dict: dict[str, Path] = {
        field: field_dir / "cutouts" for field, field_dir in zip(fields, field_dirs)
    }

    image = get_yanda(version=yanda)

    beams_col, island_col, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )

    output_dir = output_dir.absolute()
    inter_dir = output_dir / merge_name
    inter_dir.mkdir(exist_ok=True)
    data_dir = inter_dir / "cutouts"
    data_dir.mkdir(exist_ok=True)

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


def merge_parser(parent_parser: bool = False) -> argparse.ArgumentParser:
    """Merge parser.

    Args:
        parent_parser (bool, optional): Whether this is a parent parser. Defaults to False.

    Returns:
        argparse.ArgumentParser: Merge parser

    """
    # Help string to be shown using the -h option
    descStr = """
    Mosaic RACS beam fields with linmos.

    """
    # Parse the command line options
    merge_parser = argparse.ArgumentParser(
        add_help=not parent_parser,
        description=descStr,
        formatter_class=UltimateHelpFormatter,
    )
    parser = merge_parser.add_argument_group("merge arguments")
    parser.add_argument(
        "--merge_name",
        type=str,
        help="Name of the merged region",
    )

    parser.add_argument(
        "--fields",
        type=str,
        nargs="+",
        help="RACS fields to mosaic - e.g. RACS_2132-50A.",
    )

    parser.add_argument(
        "--datadirs",
        type=Path,
        nargs="+",
        help="Directories containing cutouts (in subdir outdir/cutouts)..",
    )

    parser.add_argument(
        "--output_dir",
        type=Path,
        help="Path to save merged data (in output_dir/merge_name/cutouts)",
    )
    parser.add_argument(
        "-e",
        "--epoch",
        type=int,
        default=0,
        help="Epoch of observation.",
    )

    parser.add_argument(
        "--host",
        metavar="host",
        type=str,
        default=None,
        help="Host of mongodb (probably $hostname -i).",
    )
    parser.add_argument(
        "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )
    return merge_parser


def cli():
    """Command-line interface."""
    m_parser = merge_parser(parent_parser=True)
    lin_parser = linmos_parser(parent_parser=True)

    parser = argparse.ArgumentParser(
        parents=[m_parser, lin_parser],
        formatter_class=UltimateHelpFormatter,
        description=m_parser.description,
    )
    args = parser.parse_args()

    test_db(
        host=args.host,
        username=args.username,
        password=args.password,
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
