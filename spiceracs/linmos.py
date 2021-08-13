#!/usr/bin/env python3
from pprint import pprint
from logging import disable
import os
import subprocess
import shlex
import ast
import pymongo
import numpy as np
from astropy.table import Table
from glob import glob
from astropy.coordinates import SkyCoord
import astropy.units as u
from spiceracs.utils import tqdm_dask, get_db, test_db
from IPython import embed
import astropy
import time
from dask import delayed
from dask.distributed import Client, progress, LocalCluster
from dask.diagnostics import ProgressBar
from spython.main import Client as sclient
import warnings

os.environ["OMP_NUM_THREADS"] = "1"

@delayed
def gen_seps(field, scriptdir):
    """Get beam separations

    Args:
        field (str): File name e.g. 2132-50A

    Returns:
        Table: Separation table
    """
    # get offsets
    offsets = Table.read(f"{scriptdir}/../askap_surveys/racs_low_offsets.csv")
    offsets.add_index("Beam")

    master_cat = Table.read(
        f"{scriptdir}/../askap_surveys/racs/db/epoch_0/field_data.csv"
    )
    master_cat.add_index("FIELD_NAME")
    master_cat = master_cat.loc[f"RACS_{field}"]
    if type(master_cat) is not astropy.table.row.Row:
        master_cat = master_cat[0]

    # Look for multiple SBIDs - only need one
    cats = glob(
        f"{scriptdir}/../askap_surveys/racs/db/epoch_0/beam_inf_*-RACS_{field}.csv"
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
        ra_beam = beam_coord.ra.hms
        dec_beam = beam_coord.dec.dms
        ra_field = field_coord.ra.hms
        dec_field = field_coord.dec.dms

        row = [
            beam,
            f"{offsets.loc[beam]['RA']:0.3f}",
            f"{offsets.loc[beam]['Dec']:0.3f}",
            f"{ra_beam.h:02.0f}:{ra_beam.m:02.0f}:{ra_beam.s:06.3f}",
            f"{dec_beam.d:02.0f}:{abs(dec_beam.m):02.0f}:{abs(dec_beam.s):05.2f}",
            f"{ra_field.h:02.0f}:{ra_field.m:02.0f}:{ra_field.s:06.3f}",
            f"{dec_field.d:02.0f}:{abs(dec_field.m):02.0f}:{abs(dec_field.s):05.2f}",
        ]
        cols += [row]
    tab = Table(
        data=np.array(cols), names=names, dtype=[int, str, str, str, str, str, str]
    )
    tab.add_index("BEAM")
    return tab


@delayed
def genparset(
    field,
    src_name,
    beams,
    stoke,
    datadir,
    septab,
):
    """Generate parset for LINMOS

    Args:
        field (str): Name of RACS field
        stoke (str): Stokes parameter
        datadir (str): Directory containing cutouts
        septab (Table): Table of beam seperations
        prefix (str, optional): Search for files with a prefix. Defaults to "".

    Raises:
        Exception: If no files are found in the datadir
    """
    beams = beams["beams"][field]
    ims = []
    for bm in list(set(beams['beam_list'])): # Ensure list of beams is unique!
        imfile = beams[f'{stoke.lower()}_beam{bm}_image_file']
        assert os.path.basename(os.path.dirname(
            imfile)) == src_name, "Looking in wrong directory!"
        imfile = os.path.join(os.path.abspath(datadir), imfile)
        ims.append(imfile)
    ims = sorted(ims)

    if len(ims) == 0:
        raise Exception(
            'No files found. Have you run imaging? Check your prefix?')
    imlist = "[" + ','.join([im.replace(".fits", "") for im in ims]) + "]"

    wgts = []
    for bm in list(set(beams['beam_list'])): # Ensure list of beams is unique!
        wgtsfile = beams[f'{stoke.lower()}_beam{bm}_weight_file']
        assert os.path.basename(os.path.dirname(
            wgtsfile)) == src_name, "Looking in wrong directory!"
        wgtsfile = os.path.join(os.path.abspath(datadir), wgtsfile)
        wgts.append(wgtsfile)
    wgts = sorted(wgts)

    assert len(ims) == len(wgts), "Unequal number of weights and images"

    for im, wt in zip(ims, wgts):
        assert os.path.basename(os.path.dirname(im)) == os.path.basename(
            os.path.dirname(wt)
        ), "Image and weight are in different areas!"

    weightlist = "[" + ','.join([wgt.replace(".fits", "") for wgt in wgts]) + "]"

    parset_dir = os.path.join(
        os.path.abspath(datadir), os.path.basename(os.path.dirname(ims[0]))
    )

    parset_file = os.path.join(parset_dir, f"linmos_{stoke}.in")
    parset = f"""linmos.names            = {imlist}
linmos.weights          = {weightlist}
linmos.imagetype        = fits
linmos.outname          = {ims[0][:ims[0].find('beam')]}linmos
linmos.outweight        = {wgts[0][:wgts[0].find('beam')]}linmos
linmos.weighttype       = Combined
linmos.weightstate      = Inherent
# Reference image for offsets
linmos.feeds.centre     = [{septab['FOOTPRINT_RA'][0]}, {septab['FOOTPRINT_DEC'][0]}]
linmos.feeds.spacing    = 1deg
# Beam offsets
"""
    for im in ims:
        basename = im.replace(".fits", "")
        idx = basename.find("beam")
        beamno = int(basename[len("beam") + idx : len("beam") + idx + 2])
        offset = f"linmos.feeds.{basename} = [{septab[beamno]['DELTA_RA']},{septab[beamno]['DELTA_DEC']}]\n"
        parset += offset
    with open(parset_file, "w") as f:
        f.write(parset)

    return parset_file


@delayed
def linmos(parset, fieldname, image, verbose=False):
    """Run LINMOS

    Args:
        parset (str): Parset file to run
        fieldname (str): RACS field name
        host (str): Mongo host
        verbose (bool, optional): Verbose output. Defaults to False.
    """
    workdir = os.path.dirname(parset)
    parset_name = os.path.basename(parset)
    source = os.path.basename(workdir)
    stoke = parset_name[parset_name.find(".in") - 1]
    log = parset.replace(".in", ".log")
    # os.environ["OMP_NUM_THREADS"] = "1"
    linmos_command = shlex.split(f"linmos -c {parset}")
    output = sclient.execute(image=image, command=linmos_command)
    with open(log, "w") as f:
        f.write("\n".join(output))

    new_file = glob(f"{workdir}/*.cutout.image.restored.{stoke.lower()}*.linmos.fits")

    if len(new_file) != 1:
        raise Exception(f"LINMOS file not found! -- check {log}?")

    new_file = os.path.abspath(new_file[0])
    outer = os.path.basename(os.path.dirname(new_file))
    inner = os.path.basename(new_file)
    new_file = os.path.join(outer, inner)

    if verbose:
        print(f"Cube now in {new_file}")

    query = {"Source_ID": source}
    newvalues = {"$set": {f"beams.{fieldname}.{stoke.lower()}_file": new_file}}

    return pymongo.UpdateOne(query, newvalues)

def main(
    field,
    datadir,
    client,
    host,
    username=None,
    password=None,
    dryrun=False,
    prefix="",
    stokeslist=None,
    verbose=True,
):
    """Main script
    """
    # Setup singularity image
    sclient.load('docker://csirocass/yandasoft:1.2.2-galaxy')
    image = sclient.pull(pull_folder='/tmp')

    # Use ASKAPcli to get beam separations for PB correction
    scriptdir = os.path.dirname(os.path.realpath(__file__))

    beamseps = gen_seps(field, scriptdir)

    if stokeslist is None:
        stokeslist = ["I", "Q", "U", "V"]

    if datadir is not None:
        datadir = os.path.abspath(datadir)

    cutdir = os.path.abspath(os.path.join(datadir, "cutouts"))

    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )

    # Query the DB
    query = {
        "$and": [{f"beams.{field}": {"$exists": True}}, {f"beams.{field}.DR1": True}]
    }

    island_ids = sorted(beams_col.distinct("Source_ID", query))
    big_beams = list(beams_col.find({"Source_ID": {'$in': island_ids}}).sort("Source_ID"))
    # files = sorted([name for name in glob(f"{cutdir}/*") if os.path.isdir(os.path.join(cutdir, name))])
    big_comps = list(comp_col.find({"Source_ID": {'$in': island_ids}}).sort("Source_ID"))
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
        src = beams['Source_ID']
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
                )
                parfiles.append(parfile)

    results = []
    for parset in parfiles:
        results.append(
            linmos(
                parset, field, image, verbose=True
            )
        )
    futures = client.persist(results)
    # dumb solution for https://github.com/dask/distributed/issues/4831
    time.sleep(10)
    tqdm_dask(futures, desc="Runing LINMOS", disable=(not verbose), total=len(results)*2+1)

    updates = [f.compute() for f in futures]
    if verbose:
        print("Updating database...")
    db_res = beams_col.bulk_write(updates, ordered=False)
    if verbose:
        pprint(db_res.bulk_api_result)

    print("LINMOS Done!")


def cli():
    """Command-line interface
    """
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
        "-d",
        "--dryrun",
        dest="dryrun",
        action="store_true",
        help="DON'T submit jobs (just make parsets) [False].",
    )

    parser.add_argument(
        "--prefix",
        metavar="prefix",
        type=str,
        default="",
        help="Prepend prefix to file.",
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
        host=args.host,
        username=args.username,
        password=args.password,
        verbose=verbose
    )

    main(
        field=args.field,
        datadir=args.datadir,
        client=client,
        host=args.host,
        username=args.username,
        password=args.password,
        dryrun=args.dryrun,
        prefix=args.prefix,
        stokeslist=args.stokeslist,
        verbose=verbose,
    )


if __name__ == "__main__":
    cli()
