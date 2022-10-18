#!/usr/bin/env python3
"""Create the SPICE-RACS database"""
import functools
import json
import logging as log
import os
import sys
import time
from functools import partial
from glob import glob
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import psutil
import pymongo
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord, search_around_sky
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
from IPython import embed
from tqdm import tqdm, trange

from spiceracs.utils import MyEncoder, get_db, get_field_db, getdata, test_db, yes_or_no


def source2beams(ra: float, dec: float, database: Table, max_sep=1) -> Table:
    """Find RACS beams that contain a given source position

    Args:
        ra (float): RA of source in degrees.
        dec (float): DEC of source in degrees.
        database (dict): RACS database table.
        max_sep (int, optional): Maximum seperation of source to beam centre in degrees. Defaults to 1.

    Returns:
        Table: Subset of RACS databsae table containing beams that contain the source.
    """
    c1 = SkyCoord(database["RA_DEG"] * u.deg, database["DEC_DEG"] * u.deg, frame="icrs")
    c2 = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    sep = c1.separation(c2)
    beams = database[sep < max_sep * u.deg]
    return beams


def ndix_unique(x: np.ndarray) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Find the N-dimensional array of indices of the unique values in x
    From https://stackoverflow.com/questions/54734545/indices-of-unique-values-in-array

    Args:
        x (np.ndarray): Array of values.

    Returns:
        Tuple[np.ndarray, np.ndarray]:
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


def cat2beams(
    mastercat: Table, database: Table, max_sep=1, verbose=True
) -> Tuple[np.ndarray, np.ndarray, Angle]:
    """Find the separations between sources in the master catalogue and the RACS beams

    Args:
        mastercat (Table): Master catalogue table.
        database (Table): RACS database table.
        max_sep (int, optional): Maxium source separation in degrees. Defaults to 1.
        verbose (bool, optional): Verbose output. Defaults to True.

    Returns:
        Tuple[np.ndarray, np.ndarray, Angle]: Output of astropy.coordinates.search_around_sky
    """
    log.info("Getting separations from beam centres...")
    c1 = SkyCoord(database["RA_DEG"] * u.deg, database["DEC_DEG"] * u.deg, frame="icrs")

    m_ra = mastercat["RA"]
    m_dec = mastercat["Dec"]
    if not m_ra.unit:
        m_ra = m_ra * u.deg
    if not m_dec.unit:
        m_dec = m_dec * u.deg
    c2 = SkyCoord(m_ra, m_dec, frame="icrs")

    seps = search_around_sky(c1, c2, seplimit=max_sep * u.degree)
    return seps


def source_database(
    islandcat: Table,
    compcat: Table,
    host: str,
    username: str = None,
    password: str = None,
    verbose=True,
):
    """Insert sources into the database

    Following https://medium.com/analytics-vidhya/how-to-upload-a-pandas-dataframe-to-mongodb-ffa18c0953c1

    Args:
        islandcat (Table): Island catalogue table.
        compcat (Table): Component catalogue table.
        host (str): MongoDB host IP.
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo host. Defaults to None.
        verbose (bool, optional): Verbose output. Defaults to True.
    """
    # Read in main catalogues
    # Use pandas and follow
    # https://medium.com/analytics-vidhya/how-to-upload-a-pandas-dataframe-to-mongodb-ffa18c0953c1
    df_i = islandcat.to_pandas()
    if type(df_i["Source_ID"][0]) is bytes:
        log.info("Decoding strings!")
        str_df = df_i.select_dtypes([object])
        str_df = str_df.stack().str.decode("utf-8").unstack()
        for col in str_df:
            df_i[col] = str_df[col]

    source_dict_list = df_i.to_dict("records")
    log.info("Loading islands into mongo...")
    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )
    island_col.delete_many({})  # Delete previous database
    island_col.insert_many(source_dict_list)
    count = island_col.count_documents({})
    log.info("Done loading")
    log.info(f"Total documents: {count}")

    df_c = compcat.to_pandas()
    if type(df_c["Source_ID"][0]) is bytes:
        log.info("Decoding strings!")
        str_df = df_c.select_dtypes([object])
        str_df = str_df.stack().str.decode("utf-8").unstack()
        for col in str_df:
            df_c[col] = str_df[col]

    source_dict_list = df_c.to_dict("records")

    log.info("Loading components into mongo...")
    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )
    comp_col.delete_many({})  # Delete previous database
    comp_col.insert_many(source_dict_list)
    count = comp_col.count_documents({})
    log.info("Done loading")
    log.info(f"Total documents: {count}")


def beam_database(islandcat, host, username=None, password=None, verbose=True):
    # Get pointing info from RACS database
    racs_fields = get_catalogue(verbose=verbose)

    # Get beams
    beam_list = get_beams(islandcat, racs_fields, verbose=verbose)
    log.info("Loading into mongo...")
    json_data = json.loads(json.dumps(beam_list, cls=MyEncoder))
    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )
    beams_col.delete_many({})  # Delete previous databas
    beams_col.insert_many(json_data)
    count = beams_col.count_documents({})
    log.info("Done loading")
    log.info(f"Total documents: {count}")


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
                log.warning(f"{SBID} failed...")
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


def field_database(host, username, password, verbose=True):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    basedir = f"{scriptdir}/../askap_surveys/racs/db/epoch_0"
    data_file = os.path.join(basedir, "field_data.csv")
    database = Table.read(data_file)
    df = database.to_pandas()
    field_list_dict = df.to_dict("records")
    log.info("Loading fields into mongo...")
    field_col = get_field_db(host, username=username, password=password)
    field_col.delete_many({})
    field_col.insert_many(field_list_dict)
    count = field_col.count_documents({})
    log.info("Done loading")
    log.info(f"Total documents: {count}")


def main(args, verbose=True):
    """Main script

    Arguments:
        args -- commandline args
    """

    if args.load:
        # Get database from master cat
        if args.islandcat is None:
            log.critical("Island catalogue is required!")
            islandcat = input("Enter catalogue file:")
        else:
            islandcat = args.islandcat
        if args.compcat is None:
            log.critical("Component catalogue is required!")
            compcat = input("Enter catalogue file:")
        else:
            compcat = args.compcat

        # Get the master cat
        log.info(f"Reading {islandcat}")
        island_cat = Table.read(islandcat)
        log.info(f"Reading {compcat}")
        comp_cat = Table.read(compcat)
        log.critical("This will overwrite the source database!")
        check_source = yes_or_no("Are you sure you wish to proceed?")
        log.critical("This will overwrite the beams database!")
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
                username=args.username,
                password=args.password,
                verbose=verbose,
            )
    if args.field:
        log.critical("This will overwrite the field database!")
        check_field = yes_or_no("Are you sure you wish to proceed?")
        if check_field:
            field_database(
                host=args.host,
                username=args.username,
                password=args.password,
                verbose=verbose,
            )

    else:
        log.info("Nothing to do!")


def cli():
    """Command-line interface"""
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

    parser.add_argument(
        "--field", action="store_true", help="Load field table into database [False].",
    )

    args = parser.parse_args()

    verbose = args.verbose
    if verbose:
        log.basicConfig(
            level=log.INFO,
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            force=True,
        )
    else:
        log.basicConfig(
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            force=True,
        )
    test_db(
        host=args.host, username=args.username, password=args.password, verbose=verbose
    )

    main(args, verbose=verbose)


if __name__ == "__main__":
    cli()
