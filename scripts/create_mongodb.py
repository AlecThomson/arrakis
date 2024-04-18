#!/usr/bin/env python
import argparse
import logging
import subprocess as sp
from pathlib import Path
from pprint import pformat
from typing import Optional

import pymongo
from pymongo.database import Database

from arrakis.logger import logger

logger.setLevel(logging.INFO)


class MongodError(Exception):
    pass


def start_mongod(
    dbpath: Path,
    logpath: Path,
    host: str = "localhost",
    port: Optional[int] = None,
    auth: bool = False,
):
    cmd = f"mongod --fork --dbpath {dbpath} --logpath {logpath} --bind_ip {host}"
    if auth:
        cmd += " --auth"
    if port is not None:
        cmd += f"--port {port}"
    logger.info(f"Running command: {cmd}")
    try:
        proc = sp.check_output(cmd.split())
    except sp.CalledProcessError as e:
        if e.returncode == 48:
            logger.error(
                f"mongod already running - try shutting down first with `mongod --dbpath {dbpath} --shutdown`"
            )
        logger.error(f"{e}")
        raise MongodError(f"Failed to start mongod. Command was: {cmd}")
    logger.info(proc.decode())
    logger.info("Started mongod")


def stop_mongod(
    dbpath: Path,
):
    logger.info(f"Stopping mongod with dbpath {dbpath}")
    cmd = f"mongod --dbpath {dbpath} --shutdown"
    logger.info(f"Running command: {cmd}")
    try:
        proc = sp.check_output(cmd.split())
    except sp.CalledProcessError as e:
        logger.error(f"{e}")
        raise MongodError(f"Failed to stop mongod. Command was: {cmd}")
    logger.info(proc.decode())
    logger.info("Stopped mongod")


def create_or_update_user(
    db: Database,
    username: str,
    password: str,
    roles: list,
):
    try:
        return db.command(
            "createUser",
            username,
            pwd=password,
            roles=roles,
        )
    except pymongo.errors.OperationFailure as e:
        if e.code == 51003:
            logger.warning(f"User {username} already exists - updating...")
            return db.command(
                "updateUser",
                username,
                pwd=password,
                roles=roles,
            )
        raise e


def create_admin_user(
    host: str,
    password: str,
    port: Optional[int] = None,
    username: str = "admin",
):
    logger.info(f"Creating admin user {username} on {host}:{port}")
    client = pymongo.MongoClient(host, port)
    db = client.admin
    res = create_or_update_user(
        db=db,
        username=username,
        password=password,
        roles=["userAdminAnyDatabase", "readWriteAnyDatabase"],
    )
    logger.info(pformat(res))


def create_read_only_user(
    host: str,
    password: str,
    port: Optional[int] = None,
    username: str = "reader",
):
    logger.info(f"Creating read-only user {username} on {host}:{port}")
    client = pymongo.MongoClient(host, port)
    db = client.admin
    res = create_or_update_user(
        db=db,
        username=username,
        password=password,
        roles=["readAnyDatabase"],
    )
    logger.info(pformat(res))


def main(
    dbpath: Path,
    admin_password: str,
    reader_password: str,
    host: str = "localhost",
    port: Optional[int] = None,
    admin_username: Optional[str] = "admin",
    reader_username: Optional[str] = "reader",
):
    logpath = dbpath.parent / "mongod.log"
    start_mongod(
        dbpath=dbpath,
        logpath=logpath,
        host=host,
        port=port,
        auth=False,
    )
    create_admin_user(
        host=host,
        port=port,
        password=admin_password,
        username=admin_username,
    )
    create_read_only_user(
        host=host,
        port=port,
        password=reader_password,
        username=reader_username,
    )
    stop_mongod(dbpath=dbpath)
    start_mongod(
        dbpath=dbpath,
        port=port,
        logpath=logpath,
        host=host,
        auth=True,
    )
    logger.info("Done!")


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dbpath", type=Path, required=True)
    parser.add_argument("--admin-password", type=str, required=True)
    parser.add_argument("--reader-password", type=str, required=True)
    parser.add_argument("--host", type=str, default="localhost")
    parser.add_argument("--port", type=int, default=None)
    parser.add_argument("--admin-username", type=str, default="admin")
    parser.add_argument("--reader-username", type=str, default="reader")
    args = parser.parse_args()

    main(
        dbpath=args.dbpath,
        admin_password=args.admin_password,
        reader_password=args.reader_password,
        host=args.host,
        port=args.port,
        admin_username=args.admin_username,
        reader_username=args.reader_username,
    )


if __name__ == "__main__":
    cli()
