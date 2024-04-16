#!/usr/bin/env python
import argparse
import logging
import subprocess as sp
from pathlib import Path
from pprint import pformat
from typing import Optional

import pymongo

from arrakis.logger import logger

logger.setLevel(logging.info)


class MongodError(Exception):
    pass


def start_mongod(
    dbpath: Path,
    logpath: Path,
    host: str = "localhost",
    port: Optional[int] = None,
    auth: bool = False,
):
    cmd = f"mongod --dbpath {dbpath} --logpath {logpath} --bind_ip {host} --port {port}"
    if auth:
        cmd += " --auth"
    logger.info(f"Running command: {cmd}")
    proc = sp.check_output(cmd.split(), stderr=sp.STDOUT, capture_output=True)
    logger.info(proc.stdout.decode("utf-8"))
    logger.error(proc.stderr.decode("utf-8"))
    if proc.returncode != 0:
        raise MongodError(f"Failed to start mongod. Command was: {cmd}")

    logger.info("Started mongod")


def stop_mongod(
    dbpath: Path,
):
    logger.info(f"Stopping mongod with dbpath {dbpath}")
    cmd = f"mongod --dbpath {dbpath} --shutdown"
    logger.info(f"Running command: {cmd}")
    proc = sp.check_output(cmd.split(), stderr=sp.STDOUT, capture_output=True)
    logger.info(proc.stdout.decode("utf-8"))
    logger.error(proc.stderr.decode("utf-8"))
    if proc.returncode != 0:
        raise MongodError(f"Failed to stop mongod. Command was: {cmd}")


def create_admin_user(
    host: str,
    password: str,
    port: Optional[int] = None,
    username: str = "admin",
):
    logger.info(f"Creating admin user {username} on {host}:{port}")
    client = pymongo.MongoClient(host, port)
    db = client.admin
    res = db.command(
        "createUser",
        username,
        pwd=password,
        roles=[{"role": "userAdminAnyDatabase", "db": "admin"}],
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
    res = db.command(
        "createUser",
        username,
        pwd=password,
        roles=[{"role": "userAdminAnyDatabase", "db": "admin"}],
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
