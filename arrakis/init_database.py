#!/usr/bin/env python3
"""Create the Arrakis database"""
import json
import logging
import os
import time
from pathlib import Path
from typing import Dict, List, Tuple, Union

import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord, search_around_sky
from astropy.table import Table, vstack
from pymongo.results import InsertManyResult
from tqdm import tqdm

from arrakis.logger import logger
from arrakis.utils import MyEncoder, get_db, get_field_db, test_db, yes_or_no


def source2beams(ra: float, dec: float, database: Table, max_sep: float = 1) -> Table:
    """Find RACS beams that contain a given source position

    Args:
        ra (float): RA of source in degrees.
        dec (float): DEC of source in degrees.
        database (dict): RACS database table.
        max_sep (float, optional): Maximum seperation of source to beam centre in degrees. Defaults to 1.

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
            - Array of arrays. Each array contains the indices where a given value in x is found
    """
    x_flat = x.ravel()
    ix_flat = np.argsort(x_flat)
    u, ix_u = np.unique(x_flat[ix_flat], return_index=True)
    ix_ndim = np.unravel_index(ix_flat, x.shape)
    ix_ndim = np.c_[ix_ndim] if x.ndim > 1 else ix_flat
    return u, np.split(ix_ndim, ix_u[1:])


def cat2beams(
    mastercat: Table, database: Table, max_sep: float = 1
) -> Tuple[np.ndarray, np.ndarray, Angle]:
    """Find the separations between sources in the master catalogue and the RACS beams

    Args:
        mastercat (Table): Master catalogue table.
        database (Table): RACS database table.
        max_sep (float, optional): Maxium source separation in degrees. Defaults to 1.

    Returns:
        Tuple[np.ndarray, np.ndarray, Angle]: Output of astropy.coordinates.search_around_sky
    """
    logger.info("Getting separations from beam centres...")
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
    username: Union[str, None] = None,
    password: Union[str, None] = None,
) -> Tuple[InsertManyResult, InsertManyResult]:
    """Insert sources into the database

    Following https://medium.com/analytics-vidhya/how-to-upload-a-pandas-dataframe-to-mongodb-ffa18c0953c1

    Args:
        islandcat (Table): Island catalogue table.
        compcat (Table): Component catalogue table.
        host (str): MongoDB host IP.
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo host. Defaults to None.

    Returns:
        Tuple[InsertManyResult, InsertManyResult]: Results for the islands and components inserts.
    """
    # Read in main catalogues
    # Use pandas and follow
    # https://medium.com/analytics-vidhya/how-to-upload-a-pandas-dataframe-to-mongodb-ffa18c0953c1
    df_i = islandcat.to_pandas()
    if type(df_i["Source_ID"][0]) is bytes:
        logger.info("Decoding strings!")
        str_df = df_i.select_dtypes([object])
        str_df = str_df.stack().str.decode("utf-8").unstack()
        for col in str_df:
            df_i[col] = str_df[col]

    source_dict_list = df_i.to_dict("records")
    logger.info("Loading islands into mongo...")
    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )
    island_delete_res = island_col.delete_many({})  # Delete previous database
    logger.warning(
        f"Deleted {island_delete_res.deleted_count} documents from island collection"
    )
    island_insert_res = island_col.insert_many(source_dict_list)

    count = island_col.count_documents({})
    logger.info("Done loading")
    logger.info(f"Total documents: {count}")

    df_c = compcat.to_pandas()
    if type(df_c["Source_ID"][0]) is bytes:
        logger.info("Decoding strings!")
        str_df = df_c.select_dtypes([object])
        str_df = str_df.stack().str.decode("utf-8").unstack()
        for col in str_df:
            df_c[col] = str_df[col]

    source_dict_list = df_c.to_dict("records")

    logger.info("Loading components into mongo...")
    comp_delete_res = comp_col.delete_many({})  # Delete previous database
    logger.warning(
        f"Deleted {comp_delete_res.deleted_count} documents from component collection"
    )
    comp_insert_res = comp_col.insert_many(source_dict_list)
    count = comp_col.count_documents({})
    logger.info("Done loading")
    logger.info(f"Total documents: {count}")

    return island_insert_res, comp_insert_res


def beam_database(
    database_path: Path,
    islandcat: Table,
    host: str,
    username: Union[str, None] = None,
    password: Union[str, None] = None,
    epoch: int = 0,
) -> InsertManyResult:
    """Insert beams into the database

    Args:
        islandcat (Table): Island catalogue table.
        host (str): MongoDB host IP.
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo host. Defaults to None.
        epoch (int, optional): RACS epoch to use. Defaults to 0.

    Returns:
        InsertManyResult: Result of the insert.
    """
    # Get pointing info from RACS database
    racs_fields = get_catalogue(
        survey_dir=database_path,
        epoch=epoch,
    )

    # Get beams
    beam_list = get_beams(islandcat, racs_fields)
    logger.info("Loading into mongo...")
    json_data = json.loads(json.dumps(beam_list, cls=MyEncoder))
    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )
    delete_res = beams_col.delete_many({})  # Delete previous databas
    logger.warning(f"Deleted {delete_res.deleted_count} documents from beam collection")
    insert_res = beams_col.insert_many(json_data)
    count = beams_col.count_documents({})
    logger.info("Done loading")
    logger.info(f"Total documents: {count}")

    return insert_res


def get_catalogue(survey_dir: Path, epoch: int = 0) -> Table:
    """Get the RACS catalogue for a given epoch

    Args:
        epoch (int, optional): Epoch number. Defaults to 0.

    Returns:
        Table: RACS catalogue table.

    """
    basedir = survey_dir / "racs" / "db" / f"epoch_{epoch}"
    beamfiles = list(basedir.glob("beam_inf*"))

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
            basename = beamfile.name
            idx = basename.find("RACS")
            FIELD = basename[idx:-4]
            SBID = basename[9 : idx - 1]
            try:
                tab.add_column(FIELD, name="FIELD_NAME", index=0)
                tab.add_column(int(SBID), name="SBID", index=0)
                racs_fields = vstack([racs_fields, tab])
            except TypeError:
                logger.warning(f"{SBID} failed...")
                continue
    return racs_fields


def get_beams(mastercat: Table, database: Table, epoch: int = 0) -> List[Dict]:
    """Get beams from the master catalogue

    Args:
        mastercat (Table): Master catalogue table.
        database (Table): RACS database table.

    Returns:
        List[Dict]: List of beam dictionaries.

    """
    # Get seperations on sky
    seps = cat2beams(mastercat, database, max_sep=1)
    vals, ixs = ndix_unique(seps[1])

    # Get DR1 fields
    points = np.unique(list(mastercat["Tile_ID"]))
    fields = np.array(
        [
            point.replace("_test4_1.05_","_") if epoch ==0 else point for point in points
        ]
    )

    # Fix for no 'test4' in cat
    in_dr1 = np.isin(
        [
            field.replace("_test4_1.05_","_") if epoch ==0 else field for field in database["FIELD_NAME"]
        ], 
        fields
    )

    beam_list = []
    for i, (val, idx) in enumerate(
        tqdm(zip(vals, ixs), total=len(vals), desc="Getting beams")
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
                "n_fields_DR1": np.sum([val["DR1"] for val in beam_dict.values()]),
                "beams": beam_dict,
            }
        )
    return beam_list


def field_database(
    survey_dir: Path,
    host: str,
    username: Union[str, None],
    password: Union[str, None],
    epoch: int = 0,
) -> InsertManyResult:
    """Reset and load the field database

    Args:
        host (str): Mongo host
        username (Union[str, None]): Mongo username
        password (Union[str, None]): Mongo password
        epoch (int, optional): RACS epoch number. Defaults to 0.

    Returns:
        InsertManyResult: Field insert object.
    """
    basedir = survey_dir / "racs" / "db" / f"epoch_{epoch}"
    data_file = basedir / "field_data.csv"
    database = Table.read(data_file)
    df = database.to_pandas()
    field_list_dict = df.to_dict("records")
    logger.info("Loading fields into mongo...")
    field_col = get_field_db(host, username=username, password=password)
    delete_res = field_col.delete_many({})
    logger.warning(f"Deleted documents: {delete_res.deleted_count}")
    insert_res = field_col.insert_many(field_list_dict)
    count = field_col.count_documents({})
    logger.info("Done loading")
    logger.info(f"Total documents: {count}")

    return insert_res


def main(
    load: bool = False,
    islandcat: Union[str, None] = None,
    compcat: Union[str, None] = None,
    database_path: Union[Path, None] = None,
    host: str = "localhost",
    username: Union[str, None] = None,
    password: Union[str, None] = None,
    field: bool = False,
    epoch: int = 0,
    force: bool = False,
) -> None:
    """Main script

    Args:
        load (bool, optional): Load the database. Defaults to False.
        islandcat (Union[str, None], optional): Island catalogue. Defaults to None.
        compcat (Union[str, None], optional): Component catalogue. Defaults to None.
        host (str, optional): Mongo host. Defaults to "localhost".
        username (Union[str, None], optional): Mongo username. Defaults to None.
        password (Union[str, None], optional): Mongo password. Defaults to None.
        field (bool, optional): Load the field database. Defaults to False.
        epoch (int, optional): RACS epoch to load. Defaults to 0.
        force (bool, optional): Force overwrite of database. Defaults to False.

    Raises:
        ValueError: If load is True and islandcat or compcat are None.

    """
    if force:
        logger.critical("This will overwrite the database! ALL data will be lost!")
        logger.critical("Sleeping for 30 seconds in case you want to cancel...")
        time.sleep(30)
        logger.critical("Continuing...you have been warned!")

    if load:
        # Get database from master cat
        if islandcat is None:
            logger.critical("Island catalogue is required!")
            islandcat = input("Enter catalogue file:")
        if compcat is None:
            logger.critical("Component catalogue is required!")
            compcat = input("Enter catalogue file:")
        if database_path is None:
            logger.critical("Database path is required!")
            database_path = Path(input("Enter database path:"))

        # Get the master cat
        logger.info(f"Reading {islandcat}")
        island_cat = Table.read(islandcat)
        logger.info(f"Reading {compcat}")
        comp_cat = Table.read(compcat)
        logger.critical("This will overwrite the source database!")
        check_source = (
            yes_or_no("Are you sure you wish to proceed?") if not force else True
        )
        logger.critical("This will overwrite the beams database!")
        check_beam = (
            yes_or_no("Are you sure you wish to proceed?") if not force else True
        )
        if check_source:
            source_database(
                islandcat=island_cat,
                compcat=comp_cat,
                host=host,
                username=username,
                password=password,
            )
        if check_beam:
            beam_database(
                database_path=database_path,
                islandcat=island_cat,
                host=host,
                username=username,
                password=password,
                epoch=epoch,
            )
    if field:
        logger.critical("This will overwrite the field database!")
        check_field = (
            yes_or_no("Are you sure you wish to proceed?") if not force else True
        )
        if check_field:
            field_res = field_database(
                survey_dir=database_path,
                host=host,
                username=username,
                password=password,
            )

    else:
        logger.info("Nothing to do!")

    logger.info("Done!")


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
    Arrakis Initialisation:

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
        "-u", "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "-p", "--password", type=str, default=None, help="Password of mongodb."
    )
    parser.add_argument(
        "-d", "--database-path",
        type=str,
        default=None,
        help="Path to RACS database (i.e. 'askap_surveys' repo).",
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
        "-f",
        "--field",
        action="store_true",
        help="Load field table into database [False].",
    )

    parser.add_argument(
        "-e",
        "--epoch",
        type=int,
        default=0,
        help="RACS epoch to load [0].",
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.INFO)

    test_db(host=args.host, username=args.username, password=args.password)

    main(
        database_path=Path(args.database_path) if args.database_path else None,
        load=args.load,
        islandcat=args.islandcat,
        compcat=args.compcat,
        host=args.host,
        username=args.username,
        password=args.password,
        field=args.field,
        epoch=args.epoch,
    )


if __name__ == "__main__":
    cli()
