#!/usr/bin/env python3
"""Create the Arrakis database"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord, search_around_sky
from astropy.table import Table, vstack
from pymongo.results import InsertManyResult
from tqdm import tqdm

from arrakis.logger import TqdmToLogger, UltimateHelpFormatter, logger
from arrakis.utils.database import get_beam_inf_db, get_db, get_field_db, test_db
from arrakis.utils.json import MyEncoder
from arrakis.utils.meta import yes_or_no
from arrakis.utils.pipeline import logo_str

logger.setLevel(logging.INFO)

TQDM_OUT = TqdmToLogger(logger, level=logging.INFO)


def source2beams(ra: float, dec: float, database: Table, max_sep: float = 1) -> Table:
    """Find RACS beams that contain a given source position

    Args:
        ra (float): RA of source in degrees.
        dec (float): DEC of source in degrees.
        database (dict): RACS database table.
        max_sep (float, optional): Maximum separation of source to beam centre in degrees. Defaults to 1.

    Returns:
        Table: Subset of RACS database table containing beams that contain the source.
    """
    c1 = SkyCoord(database["RA_DEG"] * u.deg, database["DEC_DEG"] * u.deg, frame="icrs")
    c2 = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    sep = c1.separation(c2)
    beams = database[sep < max_sep * u.deg]
    return beams


def ndix_unique(x: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
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
) -> tuple[np.ndarray, np.ndarray, Angle]:
    """Find the separations between sources in the master catalogue and the RACS beams

    Args:
        mastercat (Table): Master catalogue table.
        database (Table): RACS database table.
        max_sep (float, optional): Maximum source separation in degrees. Defaults to 1.

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
    epoch: int,
    username: str | None = None,
    password: str | None = None,
) -> tuple[InsertManyResult, InsertManyResult]:
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
    if isinstance(df_i["Source_ID"][0], bytes):
        logger.info("Decoding strings!")
        str_df = df_i.select_dtypes([object])
        str_df = str_df.stack().str.decode("utf-8").unstack()
        for col in str_df:
            df_i[col] = str_df[col]

    source_dict_list = df_i.to_dict("records")
    logger.info("Loading islands into mongo...")
    beams_col, island_col, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )
    island_delete_res = island_col.delete_many({})  # Delete previous database
    logger.warning(
        f"Deleted {island_delete_res.deleted_count} documents from island collection"
    )
    island_insert_res = island_col.insert_many(source_dict_list)

    count = island_col.count_documents({})
    logger.info("Done loading")
    logger.info(f"Total documents: {count}")

    logger.info("Creating index...")
    idx_res = island_col.create_index("Source_ID")
    logger.info(f"Index created: {idx_res}")

    df_c = compcat.to_pandas()
    if isinstance(df_c["Source_ID"][0], bytes):
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

    logger.info("Creating index...")
    idx_res = comp_col.create_index("Gaussian_ID")
    logger.info(f"Index created: {idx_res}")

    return island_insert_res, comp_insert_res


def beam_database(
    database_path: Path,
    islandcat: Table,
    host: str,
    epoch: int,
    username: str | None = None,
    password: str | None = None,
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
    beam_list = get_beams(islandcat, racs_fields, epoch=epoch)
    logger.info("Loading into mongo...")
    json_data = json.loads(json.dumps(beam_list, cls=MyEncoder))
    beams_col, island_col, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )
    delete_res = beams_col.delete_many({})  # Delete previous database
    logger.warning(f"Deleted {delete_res.deleted_count} documents from beam collection")
    insert_res = beams_col.insert_many(json_data)
    count = beams_col.count_documents({})
    logger.info("Done loading")
    logger.info(f"Total documents: {count}")
    logger.info("Creating index...")
    idx_res = beams_col.create_index("Source_ID")
    logger.info(f"Index created: {idx_res}")

    return insert_res


def get_catalogue(survey_dir: Path, epoch: int = 0) -> Table:
    """Get the RACS catalogue for a given epoch

    Args:
        epoch (int, optional): Epoch number. Defaults to 0.

    Returns:
        Table: RACS catalogue table.

    """
    database = read_racs_database(survey_dir, epoch, table="field_data")
    # Remove rows with SBID < 0
    database = database[database["SBID"] >= 0]

    # Init first field
    row = database[0]
    FIELD = row["FIELD_NAME"]
    SBID = row["SBID"]
    # Find FIELD and SBID in beamfile name
    racs_fields = read_racs_database(
        survey_dir, epoch, table=f"beam_inf_{SBID}-{FIELD}"
    )
    racs_fields.add_column(FIELD, name="FIELD_NAME", index=0)
    racs_fields.add_column(SBID, name="SBID", index=0)

    # Add in all others
    for row in tqdm(database[1:], desc="Reading RACS database", file=TQDM_OUT):
        try:
            tab = read_racs_database(
                survey_dir, epoch, table=f"beam_inf_{row['SBID']}-{row['FIELD_NAME']}"
            )
        except Exception as e:
            logger.error(e)
            continue
        try:
            tab.add_column(row["FIELD_NAME"], name="FIELD_NAME", index=0)
            tab.add_column(row["SBID"], name="SBID", index=0)
            racs_fields = vstack([racs_fields, tab])
        except TypeError:
            logger.error(f"{SBID} failed...")
            continue
    return racs_fields


def get_beams(mastercat: Table, database: Table, epoch: int = 0) -> list[dict]:
    """Get beams from the master catalogue

    Args:
        mastercat (Table): Master catalogue table.
        database (Table): RACS database table.

    Returns:
        List[Dict]: List of beam dictionaries.

    """
    # Get separations on sky
    seps = cat2beams(mastercat, database, max_sep=1)
    vals, ixs = ndix_unique(seps[1])

    # Get DR1 fields
    points = np.unique(list(mastercat["Tile_ID"]))
    fields = np.array(
        [
            point.replace("_test4_1.05_", "_") if epoch == 0 else point
            for point in points
        ]
    )

    # Fix for no 'test4' in cat
    in_dr1 = np.isin(
        [
            field.replace("_test4_1.05_", "_") if epoch == 0 else field
            for field in database["FIELD_NAME"]
        ],
        fields,
    )

    beam_list = []
    for i, (val, idx) in enumerate(
        tqdm(zip(vals, ixs), total=len(vals), desc="Getting beams", file=TQDM_OUT)
    ):
        beam_dict = {}
        name = mastercat[val]["Source_Name"]
        isl_id = mastercat[val]["Source_ID"]
        beams = database[seps[0][idx.astype(int)]]
        for j, field in enumerate(np.unique(beams["FIELD_NAME"])):
            ndx = beams["FIELD_NAME"] == field
            field = field.replace("_test4_1.05_", "_") if epoch == 0 else field
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


def beam_inf(
    database: Table,
    survey_dir: Path,
    host: str,
    epoch: int,
    username: str | None = None,
    password: str | None = None,
) -> InsertManyResult:
    """Get the beam information"""
    tabs: list[Table] = []
    for row in tqdm(database, desc="Reading beam info", file=TQDM_OUT):
        try:
            tab = read_racs_database(
                survey_dir=survey_dir,
                epoch=epoch,
                table=f"beam_inf_{row['SBID']}-{row['FIELD_NAME']}",
            )
        except Exception as e:
            logger.error(e)
            continue
        if len(tab) == 0:
            logger.error(f"{row['SBID']}-{row['FIELD_NAME']} failed...")
            continue
        tab.add_column(row["FIELD_NAME"], name="FIELD_NAME", index=0)
        tabs.append(tab)

    big = vstack(tabs)
    big_dict = big.to_pandas().to_dict("records")

    logger.info("Loading beam inf into mongo...")
    beam_inf_col = get_beam_inf_db(
        host, epoch=epoch, username=username, password=password
    )
    delete_res = beam_inf_col.delete_many({})
    logger.warning(f"Deleted documents: {delete_res.deleted_count}")
    insert_res = beam_inf_col.insert_many(big_dict)
    count = beam_inf_col.count_documents({})
    logger.info("Done loading")
    logger.info(f"Total documents: {count}")
    logger.info("Creating index...")
    idx_res = beam_inf_col.create_index("FIELD_NAME")
    logger.info(f"Index created: {idx_res}")

    return insert_res


def read_racs_database(
    survey_dir: Path,
    epoch: int,
    table: str,
) -> Table:
    """Read the RACS database from CSVs or postgresql

    Args:
        survey_dir (Path): Path to RACS database (i.e. 'askap_surveys/racs' repo).
        epoch (int): RACS epoch number.
        table (str): Table name.

    Returns:
        Table: RACS database table.
    """
    epoch_name = f"epoch_{epoch}"
    if survey_dir.parent.name == "postgresql:":
        logger.info("Reading RACS data from postgresql...")
        _dbstring = f"{survey_dir.parent.name}//{survey_dir.name}/{epoch_name}"
        _df = pd.read_sql(
            f'SELECT * from "{table}"',
            f"{_dbstring}",
        )
        return Table.from_pandas(_df)

    logger.info("Reading RACS data from CSVs...")
    basedir = survey_dir / "db" / epoch_name
    data_file = basedir / f"{table}.csv"
    if not data_file.exists():
        raise FileNotFoundError(f"{data_file} not found!")

    return Table.read(data_file)


def field_database(
    survey_dir: Path,
    host: str,
    epoch: int,
    username: str | None = None,
    password: str | None = None,
) -> tuple[InsertManyResult, InsertManyResult]:
    """Reset and load the field database

    Args:
        survey_dir (Path): Path to RACS database (i.e. 'askap_surveys/racs' repo).
        host (str): Mongo host
        epoch (int, optional): RACS epoch number.
        username (Union[str, None]): Mongo username
        password (Union[str, None]): Mongo password

    Returns:
        Tuple[InsertManyResult, InsertManyResult]: Field and beam info insert object.
    """
    database = read_racs_database(survey_dir, epoch, table="field_data")
    if "COMMENT" in database.colnames:
        database["COMMENT"] = database["COMMENT"].astype(str)
    # Remove rows with SBID < 0
    database = database[database["SBID"] >= 0]
    df = database.to_pandas()
    field_list_dict = df.to_dict("records")
    logger.info("Loading fields into mongo...")
    field_col = get_field_db(
        host=host, epoch=epoch, username=username, password=password
    )
    delete_res = field_col.delete_many({})
    logger.warning(f"Deleted documents: {delete_res.deleted_count}")
    insert_res = field_col.insert_many(field_list_dict)
    count = field_col.count_documents({})
    logger.info("Done loading")
    logger.info(f"Total documents: {count}")
    logger.info("Creating index...")
    idx_res = field_col.create_index("FIELD_NAME")
    logger.info(f"Index created: {idx_res}")

    beam_res = beam_inf(
        database=database,
        survey_dir=survey_dir,
        host=host,
        epoch=epoch,
        username=username,
        password=password,
    )

    return insert_res, beam_res


def main(
    load: bool = False,
    islandcat: str | None = None,
    compcat: str | None = None,
    database_path: Path | None = None,
    host: str = "localhost",
    username: str | None = None,
    password: str | None = None,
    field: bool = False,
    epochs: list[int] = 0,
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
        epochs (List[int], optional): Epochs to load. Defaults to [0].
        force (bool, optional): Force overwrite of database. Defaults to False.

    Raises:
        ValueError: If load is True and islandcat or compcat are None.

    """
    if force:
        logger.critical("This will overwrite the database! ALL data will be lost!")
        logger.critical("Sleeping for 30 seconds in case you want to cancel...")
        time.sleep(30)
        logger.critical("Continuing...you have been warned!")

    # Do checks up front
    if load:
        logger.critical("This will overwrite the source database!")
        check_source = (
            yes_or_no("Are you sure you wish to proceed?") if not force else True
        )
        logger.critical("This will overwrite the beams database!")
        check_beam = (
            yes_or_no("Are you sure you wish to proceed?") if not force else True
        )

    if field:
        logger.critical("This will overwrite the field and beam info database!")
        check_field = (
            yes_or_no("Are you sure you wish to proceed?") if not force else True
        )

    for epoch in epochs:
        logger.warning(f"Loading epoch {epoch}")
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
            logger.info(f"Reading {islandcat}. If VOTable, this may take a while...")
            island_cat = Table.read(islandcat)
            logger.info(f"Reading {compcat}. If VOTable, this may take a while...")
            comp_cat = Table.read(compcat)
            if check_source:
                source_database(
                    islandcat=island_cat,
                    compcat=comp_cat,
                    host=host,
                    epoch=epoch,
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
            if database_path is None:
                logger.critical("Database path is required!")
                database_path = Path(input("Enter database path:"))

            if check_field:
                field_res, beam_res = field_database(
                    survey_dir=database_path,
                    host=host,
                    username=username,
                    password=password,
                    epoch=epoch,
                )

        else:
            logger.info("Nothing to do!")

    logger.info("Done!")


def cli():
    """Command-line interface"""
    import argparse

    # Help string to be shown using the -h option
    descStr = f"""
    {logo_str}
    Arrakis Initialisation:

    Create MongoDB database from RACS catalogues.

    Before running make sure to start a session of mongodb e.g.
        $ mongod --dbpath=/path/to/database --bind_ip $(hostname -i)

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=UltimateHelpFormatter
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
        "-d",
        "--database-path",
        type=str,
        default=None,
        help="Path to RACS database (i.e. 'askap_surveys/racs' repo).",
    )

    parser.add_argument(
        "-i", "--islandcat", type=str, help="Master island RACS catalogue."
    )
    parser.add_argument(
        "-c", "--compcat", type=str, help="Master component RACS catalogue."
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument(
        "-l",
        "--load",
        action="store_true",
        help="Load catalogue into database.",
    )

    parser.add_argument(
        "-f",
        "--field",
        action="store_true",
        help="Load field table into database.",
    )

    parser.add_argument(
        "-e",
        "--epochs",
        type=int,
        default=[0],
        help="Epochs to load.",
        nargs="+",
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
        epochs=args.epochs,
    )


if __name__ == "__main__":
    cli()
