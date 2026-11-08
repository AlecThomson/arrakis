#!/usr/bin/env python
"""Database utilities"""

from __future__ import annotations

import warnings

import pymongo
from astropy.utils.exceptions import AstropyWarning
from pymongo.collection import Collection
from spectral_cube.utils import SpectralCubeWarning

from arrakis.logger import logger

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)


def get_sbid_observation(field_col: Collection, sbid: int) -> dict:
    """Get the field database entry for a single SBID

    Args:
        field_col (Collection): Field collection.
        sbid (int): SBID to look up.

    Raises:
        ValueError: If the SBID is not in the field database.

    Returns:
        dict: FIELD_NAME, SBID, SELECT, SCAN_START and SCAN_TINT of the observation.
    """
    field_data = field_col.find_one(
        {"SBID": sbid},
        {
            "_id": 0,
            "FIELD_NAME": 1,
            "SBID": 1,
            "SELECT": 1,
            "SCAN_START": 1,
            "SCAN_TINT": 1,
        },
    )
    if field_data is None:
        raise ValueError(f"SBID {sbid} not found in field database")
    return field_data


def get_field_observation(
    field_col: Collection, field_name: str, sbid: int | None = None
) -> dict:
    """Get the field database entry for the observation of a field

    Args:
        field_col (Collection): Field collection.
        field_name (str): Field name.
        sbid (int, optional): If given, use this SBID's observation, whatever
            its SELECT value. Otherwise use the SELECT=1 observation of the
            field. Defaults to None.

    Raises:
        ValueError: If the SBID belongs to another field, if the field has more
            than one SELECT=1 observation, or if the field is not in the database.

    Returns:
        dict: FIELD_NAME, SBID, SELECT, SCAN_START and SCAN_TINT of the observation.
    """
    # An SBID is unique, so a re-observed field's other SBIDs can't be used
    if sbid is not None:
        field_data = get_sbid_observation(field_col, sbid)
        if field_data["FIELD_NAME"] != field_name:
            raise ValueError(
                f"SBID {sbid} is field {field_data['FIELD_NAME']}, not {field_name}"
            )
        return field_data

    # SELECT=1 is the best observation of the field according to the database
    query = {"FIELD_NAME": field_name, "SELECT": 1}
    n_selected = field_col.count_documents(query)
    if n_selected > 1:
        raise ValueError(
            f"More than one SELECT=1 for {field_name} - try supplying SBID."
        )
    if n_selected == 0:
        logger.warning(f"No SELECT=1 for {field_name}, using its highest SBID.")
        query = {"FIELD_NAME": field_name}

    found = list(
        field_col.find(query, {"SBID": 1}).sort("SBID", pymongo.DESCENDING).limit(1)
    )
    if not found:
        raise ValueError(f"Field {field_name} not found in field database")
    return get_sbid_observation(field_col, found[0]["SBID"])


def validate_sbid_field_pair(field_name: str, sbid: int, field_col: Collection) -> bool:
    """Validate field and sbid pair

    Args:
        field_name (str): Field name.
        sbid (int): SBID.
        field_col (Collection): Field collection.

    Raises:
        bool: If field name and sbid pair is valid.
    """
    logger.info(f"Validating field name and SBID pair: {field_name}, {sbid}")
    field_data: dict | None = field_col.find_one({"SBID": sbid})
    if field_data is None:
        raise ValueError(f"SBID {sbid} not found in database")

    return field_data["FIELD_NAME"] == field_name


def test_db(
    host: str, username: str | None = None, password: str | None = None
) -> bool:
    """Test connection to MongoDB

    Args:
        host (str): Mongo host IP.
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo password. Defaults to None.


    Returns:
        bool: True if connection successful

    Raises:
        Exception: If connection fails.
    """
    logger.info("Testing MongoDB connection...")
    # default connection (ie, local)
    with pymongo.MongoClient(
        host=host,
        connect=False,
        username=username,
        password=password,
        authMechanism="SCRAM-SHA-256",
    ) as dbclient:  # type: pymongo.MongoClient
        try:
            dbclient.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")

        logger.info("MongoDB connection successful!")

    return True


def get_db(
    host: str,
    epoch: int,
    username: str | None = None,
    password: str | None = None,
) -> tuple[Collection, Collection, Collection]:
    """Get MongoDBs

    Args:
        host (str): Mongo host IP.
        username (str, optional): Username. Defaults to None.
        password (str, optional): Password. Defaults to None.

    Returns:
        Tuple[Collection, Collection, Collection]: beams_col, island_col, comp_col
    """
    dbclient = pymongo.MongoClient(
        host=host,
        connect=True,
        username=username,
        password=password,
        authMechanism="SCRAM-SHA-256",
    )
    mydb = dbclient[f"arrakis_epoch_{epoch}"]  # Create/open database
    comp_col = mydb["components"]  # Create/open collection
    island_col = mydb["islands"]  # Create/open collection
    beams_col = mydb["beams"]  # Create/open collection
    return beams_col, island_col, comp_col


def get_field_db(host: str, epoch: int, username=None, password=None) -> Collection:
    """Get MongoDBs

    Args:
        host (str): Mongo host IP.
        username (str, optional): Username. Defaults to None.
        password (str, optional): Password. Defaults to None.

    Returns:
        pymongo.Collection: beams_col, island_col, comp_col
    """
    dbclient = pymongo.MongoClient(
        host=host,
        connect=False,
        username=username,
        password=password,
        authMechanism="SCRAM-SHA-256",
    )  # type: pymongo.MongoClient
    mydb = dbclient[f"arrakis_epoch_{epoch}"]  # Create/open database
    field_col = mydb["fields"]  # Create/open collection
    return field_col


def get_beam_inf_db(host: str, epoch: int, username=None, password=None) -> Collection:
    """Get MongoDBs

    Args:
        host (str): Mongo host IP.
        username (str, optional): Username. Defaults to None.
        password (str, optional): Password. Defaults to None.

    Returns:
        pymongo.Collection: beams_col, island_col, comp_col
    """
    dbclient = pymongo.MongoClient(
        host=host,
        connect=False,
        username=username,
        password=password,
        authMechanism="SCRAM-SHA-256",
    )  # type: pymongo.MongoClient
    mydb = dbclient[f"arrakis_epoch_{epoch}"]  # Create/open database
    beam_inf_col = mydb["beam_inf"]  # Create/open collection
    return beam_inf_col
