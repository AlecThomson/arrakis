#!/usr/bin/env python3
"""Run LINMOS on cutouts in parallel"""
import ast
import logging
import os
import shlex
import subprocess
import sys
import time
import warnings
from glob import glob
from logging import disable
from pathlib import Path
from pprint import pformat
from typing import List, Tuple, Union

import astropy
import astropy.units as u
import dask
import numpy as np
import pkg_resources
import pymongo
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from dask import delayed, distributed
from dask.diagnostics import ProgressBar
from dask.distributed import Client, LocalCluster
from IPython import embed
from spectral_cube.utils import SpectralCubeWarning
from spython.main import Client as sclient

from arrakis.logger import logger
from arrakis.utils import chunk_dask, coord_to_string, get_db, test_db, tqdm_dask

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)

os.environ["OMP_NUM_THREADS"] = "1"


@delayed
def gen_seps(field: str, survey_dir: Path, epoch: int = 0) -> Table:
    """Get separation table for a given RACS field

    Args:
        field (str): RACS field name.

    Returns:
        Table: Table of separation for each beam.
    """
    offset_file = pkg_resources.resource_filename(
        "arrakis", f"racs_epoch_{epoch}_offsets.csv"
    )
    offsets = Table.read(offset_file)
    offsets.add_index("Beam")

    field_path = survey_dir / "db" / f"epoch_{epoch}" / "field_data.csv"
    master_cat = Table.read(field_path)
    master_cat.add_index("FIELD_NAME")
    master_cat = master_cat.loc[f"RACS_{field}"]
    if type(master_cat) is not astropy.table.row.Row:
        master_cat = master_cat[0]

    # Look for multiple SBIDs - only need one
    cats = glob(
        os.path.join(
            survey_dir, "racs", "db", "epoch_0", f"beam_inf_*-RACS_{field}.csv"
        )
    )
    beam_cat = Table.read(cats[0])
    beam_cat.add_index("BEAM_NUM")

    names = [
        "BEAM",
        "DELTA_RA",
        "DELTA_DEC",
        "BEAM_RA",
        "BEAM_DEC",
        "FOOTPRINT_RA",
        "FOOTPRINT_DEC",
    ]

    cols = []
    for beam in range(36):
        beam = int(beam)

        beam_dat = beam_cat.loc[beam]
        beam_coord = SkyCoord(beam_dat["RA_DEG"] * u.deg, beam_dat["DEC_DEG"] * u.deg)
        field_coord = SkyCoord(
            master_cat["RA_DEG"] * u.deg, master_cat["DEC_DEG"] * u.deg
        )

        beam_ra_str, beam_dec_str = coord_to_string(beam_coord)
        field_ra_str, field_dec_str = coord_to_string(field_coord)

        row = [
            beam,
            f"{offsets.loc[beam]['RA']:0.3f}",
            f"{offsets.loc[beam]['Dec']:0.3f}",
            beam_ra_str,
            beam_dec_str,
            field_ra_str,
            field_dec_str,
        ]
        cols += [row]
    tab = Table(
        data=np.array(cols), names=names, dtype=[int, str, str, str, str, str, str]
    )
    tab.add_index("BEAM")
    return tab


@delayed
def genparset(
    field: str,
    src_name: str,
    beams: dict,
    stoke: str,
    datadir: str,
    septab: Table,
    holofile: Union[str, None] = None,
) -> str:
    """Generate parset for LINMOS

    Args:
        field (str): RACS field name.
        src_name (str): RACE source name.
        beams (dict): Mongo entry for RACS beams.
        stoke (str): Stokes parameter.
        datadir (str): Data directory.
        septab (Table): Table of separations.
        holofile (str): Full path to holography file.

    Raises:
        Exception: If no files are found.

    Returns:
        str: Path to parset file.
    """
    beams = beams["beams"][field]
    ims = []
    for bm in list(set(beams["beam_list"])):  # Ensure list of beams is unique!
        imfile = beams[f"{stoke.lower()}_beam{bm}_image_file"]
        assert (
            os.path.basename(os.path.dirname(imfile)) == src_name
        ), "Looking in wrong directory!"
        imfile = os.path.join(os.path.abspath(datadir), imfile)
        ims.append(imfile)
    ims = sorted(ims)

    if len(ims) == 0:
        raise Exception("No files found. Have you run imaging? Check your prefix?")
    imlist = "[" + ",".join([im.replace(".fits", "") for im in ims]) + "]"

    wgts = []
    for bm in list(set(beams["beam_list"])):  # Ensure list of beams is unique!
        wgtsfile = beams[f"{stoke.lower()}_beam{bm}_weight_file"]
        assert (
            os.path.basename(os.path.dirname(wgtsfile)) == src_name
        ), "Looking in wrong directory!"
        wgtsfile = os.path.join(os.path.abspath(datadir), wgtsfile)
        wgts.append(wgtsfile)
    wgts = sorted(wgts)

    assert len(ims) == len(wgts), "Unequal number of weights and images"

    for im, wt in zip(ims, wgts):
        assert os.path.basename(os.path.dirname(im)) == os.path.basename(
            os.path.dirname(wt)
        ), "Image and weight are in different areas!"

    weightlist = "[" + ",".join([wgt.replace(".fits", "") for wgt in wgts]) + "]"

    parset_dir = os.path.join(
        os.path.abspath(datadir), os.path.basename(os.path.dirname(ims[0]))
    )

    parset_file = os.path.join(parset_dir, f"linmos_{stoke}.in")
    parset = f"""linmos.names            = {imlist}
linmos.weights          = {weightlist}
linmos.imagetype        = fits
linmos.outname          = {ims[0][:ims[0].find('beam')]}linmos
linmos.outweight        = {wgts[0][:wgts[0].find('beam')]}linmos
# For ASKAPsoft>1.3.0
linmos.useweightslog    = true
linmos.weighttype       = Combined
linmos.weightstate      = Inherent
"""

    if holofile is not None:
        logger.info(f"Using holography file {holofile} -- setting removeleakge to true")

        parset += f"""
linmos.primarybeam      = ASKAP_PB
linmos.primarybeam.ASKAP_PB.image = {holofile}
linmos.removeleakage    = true
"""
    else:
        logger.warning("No holography file provided - not correcting leakage!")

    with open(parset_file, "w") as f:
        f.write(parset)

    return parset_file


@delayed
def linmos(parset: str, fieldname: str, image: str, verbose=False) -> pymongo.UpdateOne:
    """Run linmos

    Args:
        parset (str): Path to parset file.
        fieldname (str): Name of RACS field.
        image (str): Name of Yandasoft image.
        verbose (bool, optional): Verbose output. Defaults to False.

    Raises:
        Exception: If LINMOS fails.
        Exception: LINMOS output not found.

    Returns:
        pymongo.UpdateOne: Mongo update object.
    """

    workdir = os.path.dirname(parset)
    rootdir = os.path.split(workdir)[0]
    junk = os.path.split(workdir)[-1]
    parset_name = os.path.basename(parset)
    source = os.path.basename(workdir)
    stoke = parset_name[parset_name.find(".in") - 1]
    log_file = parset.replace(".in", ".log")
    linmos_command = shlex.split(f"linmos -c {parset}")
    output = sclient.execute(
        image=image,
        command=linmos_command,
        bind=f"{rootdir}:{rootdir}",
        return_result=True,
    )

    outstr = "\n".join(output["message"])
    with open(log_file, "w") as f:
        f.write(outstr)
        # f.write(output['message'])

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


def main(
    field: str,
    datadir: str,
    survey_dir: Path,
    host: str,
    epoch: int = 0,
    holofile: Union[str, None] = None,
    username: Union[str, None] = None,
    password: Union[str, None] = None,
    yanda="1.3.0",
    stokeslist: Union[List[str], None] = None,
    verbose=True,
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
        stokeslist (List[str], optional): Stokes parameters to process. Defaults to None.
        verbose (bool, optional): Verbose output. Defaults to True.
    """
    # Setup singularity image
    image = get_yanda(version=yanda)

    beamseps = gen_seps(
        field=field,
        survey_dir=survey_dir,
        epoch=epoch,
    )
    if stokeslist is None:
        stokeslist = ["I", "Q", "U", "V"]

    if datadir is not None:
        datadir = os.path.abspath(datadir)

    cutdir = os.path.abspath(os.path.join(datadir, "cutouts"))

    if holofile is not None:
        holofile = os.path.abspath(holofile)

    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )
    logger.debug(f"{beams_col = }")
    # Query the DB
    query = {
        "$and": [{f"beams.{field}": {"$exists": True}}, {f"beams.{field}.DR1": True}]
    }

    island_ids = sorted(beams_col.distinct("Source_ID", query))
    big_beams = list(
        beams_col.find({"Source_ID": {"$in": island_ids}}).sort("Source_ID")
    )
    # files = sorted([name for name in glob(f"{cutdir}/*") if os.path.isdir(os.path.join(cutdir, name))])
    big_comps = list(
        comp_col.find({"Source_ID": {"$in": island_ids}}).sort("Source_ID")
    )
    comps = []
    for island_id in island_ids:
        _comps = []
        for c in big_comps:
            if c["Source_ID"] == island_id:
                _comps.append(c)
        comps.append(_comps)

    assert len(big_beams) == len(comps)

    parfiles = []
    for beams, comp in zip(big_beams, comps):
        src = beams["Source_ID"]
        if len(comp) == 0:
            warnings.warn(f"Skipping island {src} -- no components found")
            continue
        else:
            for stoke in stokeslist:
                parfile = genparset(
                    field=field,
                    src_name=src,
                    beams=beams,
                    stoke=stoke.capitalize(),
                    datadir=cutdir,
                    septab=beamseps,
                    holofile=holofile,
                )
                parfiles.append(parfile)

    results = []
    for parset in parfiles:
        results.append(linmos(parset, field, image, verbose=True))

    futures = chunk_dask(
        outputs=results,
        task_name="LINMOS",
        progress_text="Runing LINMOS",
        verbose=verbose,
    )

    updates = [f.compute() for f in futures]
    logger.info("Updating database...")
    db_res = beams_col.bulk_write(updates, ordered=False)
    logger.info(pformat(db_res.bulk_api_result))

    logger.info("LINMOS Done!")


def cli():
    """Command-line interface"""
    import argparse

    # Help string to be shown using the -h option
    descStr = f"""
    Mosaic RACS beam cubes with linmos.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
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
        "survey",
        type=str,
        help="Survey directory",
    )
    parser.add_argument(
        "--epoch",
        type=int,
        default=0,
        help="Epoch to read field data from",
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
        "-v", dest="verbose", action="store_true", help="Verbose output [False]."
    )

    parser.add_argument(
        "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )

    args = parser.parse_args()

    cluster = LocalCluster(n_workers=1)
    client = Client(cluster)

    verbose = args.verbose
    test_db(
        host=args.host, username=args.username, password=args.password, verbose=verbose
    )

    main(
        field=args.field,
        datadir=args.datadir,
        host=args.host,
        holofile=args.holofile,
        username=args.username,
        password=args.password,
        yanda=args.yanda,
        stokeslist=args.stokeslist,
        verbose=verbose,
    )


if __name__ == "__main__":
    cli()
