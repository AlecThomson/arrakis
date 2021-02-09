#!/usr/bin/env python3
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
from tqdm import tqdm, trange
from IPython import embed
import astropy
import dask
from dask import delayed
from dask.distributed import Client, progress, LocalCluster
from dask.diagnostics import ProgressBar


@delayed
def gen_seps(field, scriptdir):
    """Get beam separations

    Args:
        field (str): File name e.g. 2132-50A

    Returns:
        Table: Separation table
    """
    # get offsets
    offsets = Table.read(f'{scriptdir}/../askap_surveys/racs/offsets.csv')
    offsets.add_index('Beam')

    master_cat = Table.read(
        f'{scriptdir}/../askap_surveys/racs/db/epoch_0/field_data.csv')
    master_cat.add_index('FIELD_NAME')
    master_cat = master_cat.loc[f"RACS_{field}"]
    if type(master_cat) is not astropy.table.row.Row:
        master_cat = master_cat[0]

    # Look for multiple SBIDs - only need one
    cats = glob(
        f'{scriptdir}/../askap_surveys/racs/db/epoch_0/beam_inf_*-RACS_{field}.csv')
    beam_cat = Table.read(cats[0])
    beam_cat.add_index('BEAM_NUM')

    names = ["BEAM", "DELTA_RA", "DELTA_DEC", "BEAM_RA",
             "BEAM_DEC", "FOOTPRINT_RA", "FOOTPRINT_DEC"]

    cols = []
    for beam in range(36):
        beam = int(beam)

        beam_dat = beam_cat.loc[beam]
        beam_coord = SkyCoord(
            beam_dat['RA_DEG']*u.deg, beam_dat['DEC_DEG']*u.deg)
        field_coord = SkyCoord(
            master_cat['RA_DEG']*u.deg, master_cat['DEC_DEG']*u.deg)
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
    tab = Table(data=np.array(cols), names=names, dtype=[int, str,str,str,str,str,str])
    tab.add_index("BEAM")
    return tab


@delayed
def genparset(field, stoke, datadir, septab, prefix=""):
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
    ims = sorted(
        glob(
            f"{datadir}/*.cutout.image.restored.{stoke.lower()}.*.beam*[00-35.total.fits]"
        )
    )
    if len(ims) == 0:
        print(
            f"{datadir}/*.cutout.image.restored.{stoke.lower()}.*.beam*[00-35.total.fits]"
        )
        raise Exception(
            'No files found. Have you run imaging? Check your prefix?')
    imlist = "["
    for im in ims:
        imlist += im.replace(".fits", "") + ","
    imlist = imlist[:-1] + "]"

    wgts = sorted(
        glob(
            f"{datadir}/*.weights.{stoke.lower()}.*.beam*[00-35.fits]"
        )
    )
    weightlist = "["
    for wgt in wgts:
        weightlist += wgt.replace(".fits", "") + ","
    weightlist = weightlist[:-1] + "]"

    parset_dir = datadir

    parset_file = f"{parset_dir}/linmos_{stoke}.in"
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
        basename = im.replace('.fits', '')
        idx = basename.find('beam')
        beamno = int(basename[len('beam')+idx:len('beam')+idx+2])
        offset = f"linmos.feeds.{basename} = [{septab[beamno]['DELTA_RA']},{septab[beamno]['DELTA_DEC']}]\n"
        parset += offset
    with open(parset_file, "w") as f:
        f.write(parset)

    return parset_file


@delayed
def linmos(parset, fieldname, host, verbose=False):
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
    stoke = parset_name[parset_name.find('.in')-1]
    log = parset.replace('.in', '.log')

    # os.chdir(workdir)
    linmos_command = shlex.split(f"linmos -c {parset}")
    output = subprocess.run(linmos_command, capture_output=True)
    with open(log, 'w') as f:
        f.write(output.stdout.decode("utf-8"))


    new_file = glob(
        f'{workdir}/*.cutout.image.restored.{stoke.lower()}*.linmos.fits'
        )

    assert len(new_file) == 1, "LINMOS file not found!"

    if verbose:
        print(f'Cube now in {new_file[0]}')

    with pymongo.MongoClient(host=host, connect=False) as dbclient:
        # default connection (ie, local)
        mydb = dbclient['spiceracs']  # Create/open database
        beams_col = mydb['beams']  # Create/open collection

    query = {"Source_ID": source}
    newvalues = {
        "$set":
        {
            f"beams.{fieldname}.{stoke.lower()}_file": new_file[0]
        }
    }

    beams_col.update_one(query, newvalues)


def main(field, datadir, client, host, dryrun=False, prefix="", stokeslist=None, verbose=True):
    """Main script
    """

    # Use ASKAPcli to get beam separations for PB correction
    scriptdir = os.path.dirname(os.path.realpath(__file__))

    beamseps = gen_seps(field, scriptdir)

    if stokeslist is None:
        stokeslist = ["I", "Q", "U", "V"]

    if datadir is not None:
        if datadir[-1] == '/':
            datadir = datadir[:-1]

    cutdir = f"{datadir}/cutouts"
    files = sorted([name for name in glob(f"{cutdir}/*") if os.path.isdir(os.path.join(cutdir, name))])

    parfiles = []
    for file in files:
        if os.path.basename(file) == 'slurmFiles':
            continue
        for stoke in stokeslist:
            parfile = genparset(field, stoke.capitalize(),
                                file, beamseps, prefix=prefix)
            parfiles.append(parfile)

    results = []
    for parset in parfiles:
        results.append(linmos(parset, field, host, verbose=True))


    results = client.persist(results)
    if verbose:
            print("Running LINMOS...")
    progress(results)


    print('LINMOS Done!')


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
        description=descStr,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "field",
        metavar="field",
        type=str,
        help="RACS field to mosaic - e.g. 2132-50A."
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
        nargs='+',
        type=str,
        help="List of Stokes parameters to image [ALL]",
    )

    parser.add_argument(
        'host',
        metavar='host',
        type=str,
        help='Host of mongodb (probably $hostname -i).'
    )

    parser.add_argument(
        "-v",
        dest="verbose",
        action="store_true",
        help="Verbose output [False]."
    )

    args = parser.parse_args()

    cluster = LocalCluster(n_workers=20)
    client = Client(cluster)

    verbose = args.verbose
    if verbose:
        print('Testing MongoDB connection...')
    # default connection (ie, local)
    with pymongo.MongoClient(host=args.host, connect=False) as dbclient:
        try:
            dbclient.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")
        else:
            if verbose:
                print('MongoDB connection succesful!')

    main(field=args.field,
         datadir=args.datadir,
         client=client,
         host=args.host,
         dryrun=args.dryrun,
         prefix=args.prefix,
         stokeslist=args.stokeslist,
         verbose=verbose
         )


if __name__ == "__main__":
    cli()
