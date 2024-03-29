#!/usr/bin/env python3
"""Run LINMOS on cutouts in parallel"""
import logging
import os
import shlex
import warnings
from glob import glob
from pathlib import Path
from pprint import pformat
from typing import Dict, List
from typing import NamedTuple as Struct
from typing import Optional

import pymongo
from astropy.utils.exceptions import AstropyWarning
from prefect import flow, task, unmapped
from racs_tools import beamcon_3D
from spectral_cube.utils import SpectralCubeWarning
from spython.main import Client as sclient

from arrakis.logger import logger
from arrakis.utils.database import get_db, test_db

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)

os.environ["OMP_NUM_THREADS"] = "1"

logger.setLevel(logging.INFO)


class ImagePaths(Struct):
    """Class to hold image paths"""

    images: List[Path]
    """List of image paths"""
    weights: List[Path]
    """List of weight paths"""


@task(name="Find images")
def find_images(
    field: str,
    beams: dict,
    stoke: str,
    datadir: Path,
) -> ImagePaths:
    """Find the images and weights for a given field and stokes parameter

    Args:
        field (str): Field name.
        beams (dict): Beam information.
        stoke (str): Stokes parameter.
        datadir (Path): Data directory.

    Raises:
        Exception: If no files are found.

    Returns:
        ImagePaths: List of images and weights.
    """
    logger.setLevel(logging.INFO)

    src_name = beams["Source_ID"]
    field_beams = beams["beams"][field]

    # First check that the images exist
    image_list: List[Path] = []
    for bm in list(set(field_beams["beam_list"])):  # Ensure list of beams is unique!
        imfile = Path(field_beams[f"{stoke.lower()}_beam{bm}_image_file"])
        assert imfile.parent.name == src_name, "Looking in wrong directory!"
        new_imfile = datadir.resolve() / imfile
        image_list.append(new_imfile)
    image_list = sorted(image_list)

    if len(image_list) == 0:
        raise Exception("No files found. Have you run imaging? Check your prefix?")

    weight_list: List[Path] = []
    for bm in list(set(field_beams["beam_list"])):  # Ensure list of beams is unique!
        wgtsfile = Path(field_beams[f"{stoke.lower()}_beam{bm}_weight_file"])
        assert wgtsfile.parent.name == src_name, "Looking in wrong directory!"
        new_wgtsfile = datadir.resolve() / wgtsfile
        weight_list.append(new_wgtsfile)
    weight_list = sorted(weight_list)

    assert len(image_list) == len(weight_list), "Unequal number of weights and images"

    for im, wt in zip(image_list, weight_list):
        assert (
            im.parent.name == wt.parent.name
        ), "Image and weight are in different areas!"

    return ImagePaths(image_list, weight_list)


@task(name="Smooth images")
def smooth_images(
    image_dict: Dict[str, ImagePaths],
) -> Dict[str, ImagePaths]:
    """Smooth cubelets to a common resolution

    Args:
        image_list (ImagePaths): List of cubelets to smooth.

    Returns:
        ImagePaths: Smoothed cubelets.
    """
    smooth_dict: Dict[str, ImagePaths] = {}
    for stoke, image_list in image_dict.items():
        infiles: List[str] = []
        for im in image_list.images:
            if im.suffix == ".fits":
                infiles.append(im.resolve().as_posix())
        datadict = beamcon_3D.main(
            infile=[im.resolve().as_posix() for im in image_list.images],
            uselogs=False,
            mode="total",
            conv_mode="robust",
            suffix="cres",
        )
        smooth_files: List[Path] = []
        for key, val in datadict.items():
            smooth_files.append(Path(val["outfile"]))
        smooth_dict[stoke] = ImagePaths(smooth_files, image_list.weights)

    return smooth_dict


@task(name="Generate parset")
def genparset(
    image_paths: ImagePaths,
    stoke: str,
    datadir: Path,
    holofile: Optional[Path] = None,
) -> str:
    """Generate parset for LINMOS

    Args:
        image_paths (ImagePaths): List of images and weights.
        stoke (str): Stokes parameter.
        datadir (Path): Data directory.
        holofile (Path, optional): Path to the holography file to include in the bind list. Defaults to None.

    Raises:
        Exception: If no files are found.

    Returns:
        str: Path to parset file.
    """
    logger.setLevel(logging.INFO)

    image_string = f"[{','.join([im.resolve().with_suffix('').as_posix() for im in image_paths.images])}]"
    weight_string = f"[{','.join([im.resolve().with_suffix('').as_posix() for im in image_paths.weights])}]"

    parset_dir = datadir.resolve() / image_paths.images[0].parent.name

    first_image = image_paths.images[0].resolve().with_suffix("").as_posix()
    first_weight = image_paths.weights[0].resolve().with_suffix("").as_posix()
    linmos_image_str = f"{first_image[:first_image.find('beam')]}linmos"
    linmos_weight_str = f"{first_weight[:first_weight.find('beam')]}linmos"

    parset_file = os.path.join(parset_dir, f"linmos_{stoke}.in")
    parset = f"""linmos.names            = {image_string}
linmos.weights          = {weight_string}
linmos.imagetype        = fits
linmos.outname          = {linmos_image_str}
linmos.outweight        = {linmos_weight_str}
# For ASKAPsoft>1.3.0
linmos.useweightslog    = true
linmos.weighttype       = Combined
linmos.weightstate      = Inherent
"""

    if holofile is not None:
        logger.info(f"Using holography file {holofile} -- setting removeleakge to true")

        parset += f"""
linmos.primarybeam      = ASKAP_PB
linmos.primarybeam.ASKAP_PB.image = {holofile.resolve().as_posix()}
linmos.removeleakage    = true
"""
    else:
        logger.warning("No holography file provided - not correcting leakage!")

    with open(parset_file, "w") as f:
        f.write(parset)

    return parset_file


@task(name="Run linmos")
def linmos(
    parset: Optional[str], fieldname: str, image: str, holofile: Path
) -> Optional[pymongo.UpdateOne]:
    """Run linmos

    Args:
        parset (str): Path to parset file.
        fieldname (str): Name of RACS field.
        image (str): Name of Yandasoft image.
        holofile (Path): Path to the holography file to include in the bind list.
        verbose (bool, optional): Verbose output. Defaults to False.

    Raises:
        Exception: If LINMOS fails.
        Exception: LINMOS output not found.

    Returns:
        pymongo.UpdateOne: Mongo update object.
    """
    logger.setLevel(logging.INFO)

    if parset is None:
        return

    workdir = os.path.dirname(parset)
    rootdir = os.path.split(workdir)[0]
    parset_name = os.path.basename(parset)
    source = os.path.basename(workdir)
    stoke = parset_name[parset_name.find(".in") - 1]
    log_file = parset.replace(".in", ".log")
    linmos_command = shlex.split(f"linmos -c {parset}")

    holo_folder = holofile.parent

    output = sclient.execute(
        image=image,
        command=linmos_command,
        bind=f"{rootdir}:{rootdir},{holo_folder}:{holo_folder}",
        return_result=True,
    )

    outstr = "\n".join(output["message"])
    with open(log_file, "w") as f:
        f.write(outstr)

    if output["return_code"] != 0:
        raise Exception(f"LINMOS failed! Check '{log_file}'")

    new_files = glob(f"{workdir}/*.cutout.image.restored.{stoke.lower()}*.linmos.fits")

    if len(new_files) != 1:
        raise Exception(f"LINMOS file not found! -- check {log_file}?")

    new_file = os.path.abspath(new_files[0])
    outer = os.path.basename(os.path.dirname(new_file))
    inner = os.path.basename(new_file)
    new_file = os.path.join(outer, inner)

    logger.info(f"Cube now in {workdir}/{inner}")

    query = {"Source_ID": source}
    newvalues = {"$set": {f"beams.{fieldname}.{stoke.lower()}_file": new_file}}

    return pymongo.UpdateOne(query, newvalues)


def get_yanda(version="1.3.0") -> str:
    """Pull yandasoft image from dockerhub.

    Args:
        version (str, optional): Yandasoft version. Defaults to "1.3.0".

    Returns:
        str: Path to yandasoft image.
    """
    sclient.load(f"docker://csirocass/yandasoft:{version}-galaxy")
    image = os.path.abspath(sclient.pull())
    return image


@flow(name="LINMOS")
def main(
    field: str,
    datadir: Path,
    host: str,
    epoch: int,
    holofile: Optional[Path] = None,
    username: Optional[str] = None,
    password: Optional[str] = None,
    yanda: str = "1.3.0",
    yanda_img: Optional[Path] = None,
    stokeslist: Optional[List[str]] = None,
    limit: Optional[int] = None,
) -> None:
    """Main script

    Args:
        field (str): RACS field name.
        datadir (str): Data directory.
        host (str): MongoDB host IP.
        holofile (str): Path to primary beam file.
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo password. Defaults to None.
        yanda (str, optional): Yandasoft version. Defaults to "1.3.0".
        yanda_img (Path, optional): Path to a yandasoft singularirt image. If `None`, the container version `yanda` will be downloaded. Defaults to None.
        stokeslist (List[str], optional): Stokes parameters to process. Defaults to None.
        verbose (bool, optional): Verbose output. Defaults to True.
    """
    # Setup singularity image
    image = get_yanda(version=yanda) if yanda_img is None else yanda_img

    logger.info(f"The yandasoft image is {image=}")

    if stokeslist is None:
        stokeslist = ["I", "Q", "U", "V"]

    cutdir = datadir / "cutouts"

    beams_col, island_col, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )
    logger.debug(f"{beams_col = }")
    # Query the DB
    query = {"$and": [{f"beams.{field}": {"$exists": True}}]}

    logger.info(f"The query is {query=}")

    island_ids: List[str] = sorted(beams_col.distinct("Source_ID", query))
    big_beams: List[dict] = list(
        beams_col.find({"Source_ID": {"$in": island_ids}}).sort("Source_ID")
    )
    big_comps: List[dict] = list(
        comp_col.find({"Source_ID": {"$in": island_ids}}).sort("Source_ID")
    )
    comps: List[List[dict]] = []
    for island_id in island_ids:
        _comps = []
        for c in big_comps:
            if c["Source_ID"] == island_id:
                _comps.append(c)
        comps.append(_comps)

    assert len(big_beams) == len(comps)

    if limit is not None:
        logger.critical(f"Limiting to {limit} islands")
        big_beams = big_beams[:limit]
        comps = comps[:limit]

    all_parfiles = []
    for stoke in stokeslist:
        image_paths = find_images.map(
            field=unmapped(field),
            beams=big_beams,
            stoke=unmapped(stoke.capitalize()),
            datadir=unmapped(cutdir),
        )
        parfiles = genparset.map(
            image_paths=image_paths,
            stoke=unmapped(stoke.capitalize()),
            datadir=unmapped(cutdir),
            holofile=unmapped(holofile),
        )
        all_parfiles.extend(parfiles)

    results = linmos.map(
        all_parfiles,
        unmapped(field),
        unmapped(str(image)),
        unmapped(holofile),
    )

    updates = [f.result() for f in results]
    updates = [u for u in updates if u is not None]
    logger.info("Updating database...")
    db_res = beams_col.bulk_write(updates, ordered=False)
    logger.info(pformat(db_res.bulk_api_result))

    logger.info("LINMOS Done!")


def cli():
    """Command-line interface"""
    import argparse

    # Help string to be shown using the -h option
    descStr = """
    Mosaic RACS beam cubes with linmos.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field to mosaic - e.g. 2132-50A."
    )

    parser.add_argument(
        "datadir",
        metavar="datadir",
        type=str,
        help="Directory containing cutouts (in subdir outdir/cutouts)..",
    )
    parser.add_argument(
        "--holofile", type=str, default=None, help="Path to holography image"
    )

    parser.add_argument(
        "--yanda",
        type=str,
        default="1.3.0",
        help="Yandasoft version to pull from DockerHub [1.3.0].",
    )
    parser.add_argument(
        "--yanda_image",
        default=None,
        type=Path,
        help="Path to an existing yandasoft singularity container image. ",
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
        "-v", dest="verbose", action="store_true", help="Verbose output [False]."
    )

    parser.add_argument(
        "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )
    parser.add_argument(
        "--limit",
        type=Optional[int],
        default=None,
        help="Limit the number of islands to process.",
    )

    args = parser.parse_args()

    verbose = args.verbose
    test_db(
        host=args.host, username=args.username, password=args.password, verbose=verbose
    )

    main(
        field=args.field,
        datadir=Path(args.datadir),
        host=args.host,
        epoch=args.epoch,
        holofile=Path(args.holofile),
        username=args.username,
        password=args.password,
        yanda=args.yanda,
        yanda_img=args.yanda_image,
        stokeslist=args.stokeslist,
        limit=args.limit,
    )


if __name__ == "__main__":
    cli()
