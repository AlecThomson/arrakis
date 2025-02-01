#!/usr/bin/env python
"""I/O utilities"""

from __future__ import annotations

import logging
import os
import shlex
import subprocess as sp
import warnings
from pathlib import Path

from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from spectral_cube.utils import SpectralCubeWarning

from arrakis.logger import TqdmToLogger, logger

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)

TQDM_OUT = TqdmToLogger(logger, level=logging.INFO)


def verify_tarball(
    tarball: str | Path,
):
    cmd = f"tar -tvf {tarball}"
    logger.info(f"Verifying tarball {tarball}")
    popen = sp.Popen(shlex.split(cmd), stderr=sp.PIPE)
    with popen.stderr:
        for line in iter(popen.stderr.readline, b""):
            logger.error(line.decode().strip())
    exitcode = popen.wait()
    return exitcode == 0


def parse_env_path(env_path: Path | str) -> Path:
    """Parse an environment path.

    Args:
        env_path (str): Environment path.

    Returns:
        Path: Parsed path.
    """
    if isinstance(env_path, Path):
        env_path_str = env_path.as_posix()
    else:
        env_path_str = env_path
    return Path(os.path.expandvars(env_path_str))


def gettable(tabledir: str, keyword: str) -> tuple[Table, str]:
    """Get a table from a directory given a keyword to glob.

    Args:
        tabledir (str): Directory.
        keyword (str): Keyword to glob for.

    Returns:
        Tuple[Table, str]: Table and it's file location.
    """
    table_path = Path(tabledir)
    # Glob out the necessary files
    files = list(table_path.glob(f"*.{keyword}*.xml"))  # Selvay VOTab
    filename = files[0]
    logger.info(f"Getting table data from {filename}...")

    # Get selvay data from VOTab
    table = Table.read(filename, format="votable")
    table = table.to_pandas()
    str_df = table.select_dtypes([object])
    str_df = str_df.stack().str.decode("utf-8").unstack()
    for col in str_df:
        table[col] = str_df[col]
    return table, filename.as_posix()
