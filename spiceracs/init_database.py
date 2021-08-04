#!/usr/bin/env python3
from spiceracs.utils import get_db, test_db
from IPython import embed
from functools import partial
import functools
import psutil
import pymongo
import json
import time
from tqdm import tqdm, trange
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
import sys
import os
from glob import glob
from astropy.table import Table, vstack
from spiceracs.utils import getdata, MyEncoder


def yes_or_no(question):
    while "Please answer 'y' or 'n'":
        reply = str(input(question + " (y/n): ")).lower().strip()
        if reply[:1] == "y":
            return True
        if reply[:1] == "n":
            return False


def source2beams(ra, dec, database, max_sep=1):
    """Find RACS beams containing a position.

    Arguments:
        ra {float} -- RA of source in degrees.
        dec {float} -- DEC of source in degrees.
        database {astropy.table.table.Table} -- RACS database loaded as one Table.
    Keyword Arguments:
        max_sep {int} -- Maximum angular distance to beam centre in degrees (default: {4})

    Returns:
        beams {astropy.table.table.Table} -- RACS database rows matching the source location.
    """
    c1 = SkyCoord(database["RA_DEG"] * u.deg, database["DEC_DEG"] * u.deg, frame="icrs")
    c2 = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    sep = c1.separation(c2)
    beams = database[sep < max_sep * u.deg]
    return beams


def ndix_unique(x):
    """
    From: https://stackoverflow.com/questions/54734545/indices-of-unique-values-in-array
    Returns an N-dimensional array of indices
    of the unique values in x
    ----------
    x: np.array
       Array with arbitrary dimensions
    Returns
    -------
    - 1D-array of sorted unique values
    - Array of arrays. Each array contains the indices where a
      given value in x is found
    """
    x_flat = x.ravel()
    ix_flat = np.argsort(x_flat)
    u, ix_u = np.unique(x_flat[ix_flat], return_index=True)
    ix_ndim = np.unravel_index(ix_flat, x.shape)
    ix_ndim = np.c_[ix_ndim] if x.ndim > 1 else ix_flat
    return u, np.split(ix_ndim, ix_u[1:])


def cat2beams(mastercat, database, max_sep=1, verbose=True):
    if verbose:
        print("Getting separations from beam centres...")
    c1 = SkyCoord(database["RA_DEG"] * u.deg, database["DEC_DEG"] * u.deg, frame="icrs")
    c2 = SkyCoord(mastercat["RA"] * u.deg, mastercat["Dec"] * u.deg, frame="icrs")

    seps = search_around_sky(c1, c2, seplimit=max_sep * u.degree)
    return seps


def source_database(
    islandcat, compcat, host, username=None, password=None, verbose=True
):
    # Read in main catalogues
    # Use pandas and follow
    # https://medium.com/analytics-vidhya/how-to-upload-a-pandas-dataframe-to-mongodb-ffa18c0953c1
    islandcat = islandcat.to_pandas()
    str_df = islandcat.select_dtypes([np.object])
    str_df = str_df.stack().str.decode("utf-8").unstack()
    for col in str_df:
        islandcat[col] = str_df[col]

    source_dict_list = islandcat.to_dict("records")
    if verbose:
        print("Loading islands into mongo...")
    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )
    island_col.delete_many({})  # Delete previous database
    island_col.insert_many(source_dict_list)
    count = island_col.count_documents({})
    if verbose:
        print("Done loading")
        print("Total documents:", count)

    compcat = compcat.to_pandas()
    str_df = compcat.select_dtypes([np.object])
    str_df = str_df.stack().str.decode("utf-8").unstack()
    for col in str_df:
        compcat[col] = str_df[col]

    source_dict_list = compcat.to_dict("records")

    if verbose:
        print("Loading components into mongo...")
    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )
    comp_col.delete_many({})  # Delete previous database
    comp_col.insert_many(source_dict_list)
    count = comp_col.count_documents({})
    if verbose:
        print("Done loading")
        print("Total documents:", count)


def beam_database(islandcat, host, username=None, password=None, verbose=True):
    # Get pointing info from RACS database
    racs_fields = get_catalogue(verbose=verbose)

    # Get beams
    beam_list = get_beams(islandcat, racs_fields, verbose=verbose)
    if verbose:
        print("Loading into mongo...")
    json_data = json.loads(json.dumps(beam_list, cls=MyEncoder))
    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )
    beams_col.delete_many({})  # Delete previous databas
    beams_col.insert_many(json_data)
    count = beams_col.count_documents({})
    if verbose:
        print("Done loading")
        print("Total documents:", count)


def get_catalogue(verbose=True):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    basedir = f"{scriptdir}/../askap_surveys/racs/db/epoch_0"
    beamfiles = glob(f"{basedir}/beam_inf*")

    # Init first field
    beamfile = beamfiles[0]
    racs_fields = Table.read(beamfile)
    basename = os.path.basename(beamfile)
    idx = basename.find("RACS")
    FIELD = basename[idx:-4]
    SBID = basename[9 : idx - 1]
    racs_fields.add_column(FIELD, name="FIELD_NAME", index=0)
    racs_fields.add_column(int(SBID), name="SBID", index=0)

    # Add in all others
    for i, beamfile in enumerate(tqdm(beamfiles, desc="Reading RACS database")):
        if i == 0:
            continue
        else:
            tab = Table.read(beamfile)
            basename = os.path.basename(beamfile)
            idx = basename.find("RACS")
            FIELD = basename[idx:-4]
            SBID = basename[9 : idx - 1]
            try:
                tab.add_column(FIELD, name="FIELD_NAME", index=0)
                tab.add_column(int(SBID), name="SBID", index=0)
                racs_fields = vstack([racs_fields, tab])
            except TypeError:
                print(f"{SBID} failed...")
                continue
    return racs_fields


def get_beams(mastercat, database, verbose=True):
    # Get seperations on sky
    seps = cat2beams(mastercat, database, max_sep=1, verbose=verbose)
    vals, ixs = ndix_unique(seps[1])

    # Get DR1 fields
    points = np.unique(list(mastercat["Tile_ID"]))
    fields = np.array([point[-8:] for point in points])

    # Fix for no 'test4' in cat
    # in_dr1 = np.isin(database['FIELD_NAME'], points)
    in_dr1 = np.isin([field[-8:] for field in database["FIELD_NAME"]], fields)

    beam_list = []
    for i, (val, idx) in enumerate(
        tqdm(
            zip(vals, ixs), total=len(vals), desc="Getting beams", disable=(not verbose)
        )
    ):
        beam_dict = {}
        ra = mastercat[val]["RA"]
        dec = dec = mastercat[val]["Dec"]
        name = mastercat[val]["Source_Name"]
        isl_id = mastercat[val]["Source_ID"]
        beams = database[seps[0][idx.astype(int)]]
        for j, field in enumerate(np.unique(beams["FIELD_NAME"])):
            ndx = beams["FIELD_NAME"] == field
            field = field[-8:]
            beam_dict.update(
                {
                    field: {
                        "beam_list": list(beams["BEAM_NUM"][ndx]),
                        "SBIDs": list(np.unique(beams["SBID"][ndx])),
                        "DR1": bool(np.unique(in_dr1[seps[0][idx.astype(int)]][ndx])),
                    }
                }
            )

        beam_list.append(
            {
                "Source_Name": name,
                "Source_ID": isl_id,
                "n_fields": len(beam_dict.keys()),
                "n_fields_DR1": sum([val["DR1"] for val in beam_dict.values()]),
                "beams": beam_dict,
            }
        )
    return beam_list


def main(args, verbose=True):
    """Main script

    Arguments:
        args {[type]} -- commandline args
    """

    if args.load:
        # Get database from master cat
        if args.islandcat is None:
            print("Island catalogue is required!")
            islandcat = input("Enter catalogue file:")
        else:
            islandcat = args.islandcat
        if args.compcat is None:
            print("Component catalogue is required!")
            compcat = input("Enter catalogue file:")
        else:
            compcat = args.compcat

        # Get the master cat
        island_cat = Table.read(islandcat)
        comp_cat = Table.read(compcat)
        print("This will overwrite the source database!")
        check_source = yes_or_no("Are you sure you wish to proceed?")
        print("This will overwrite the beams database!")
        check_beam = yes_or_no("Are you sure you wish to proceed?")
        if check_source:
            source_database(
                islandcat=island_cat,
                compcat=comp_cat,
                host=args.host,
                username=args.username,
                password=args.password,
                verbose=verbose,
            )
        if check_beam:
            beam_database(
                islandcat=island_cat, 
                host=args.host,
                username=args.host,
                password=args.password, 
                verbose=verbose
                )

    else:
        print("Nothing to do!")


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
    SPICE-RACS Initialisation:
    
    Create MongoDB database from RACS catalogues.

    Before running make sure to start a session of mongodb e.g.
        $ mongod --dbpath=/path/to/database --bind_ip $(hostname -i)

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
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
        "-i", "--islandcat", type=str, help="Master island RACS catalogue."
    )
    parser.add_argument(
        "-c", "--compcat", type=str, help="Master component RACS catalogue."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose output [False]."
    )
    parser.add_argument(
        "-l",
        "--load",
        action="store_true",
        help="Load catalogue into database [False].",
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
