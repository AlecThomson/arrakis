#!/usr/bin/env python
"""I/O utilities"""

from __future__ import annotations

import logging
import os
import shlex
import stat
import subprocess as sp
import warnings
from pathlib import Path

from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from spectral_cube.utils import SpectralCubeWarning
from tqdm.auto import tqdm

from arrakis.logger import TqdmToLogger, logger
from arrakis.utils.exceptions import SameFileError, SpecialFileError
from arrakis.utils.typing import PathLike

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


def parse_env_path(env_path: PathLike) -> Path:
    """Parse an environment path.

    Args:
        env_path (str): Environment path.

    Returns:
        Path: Parsed path.
    """
    if isinstance(env_path, Path):
        env_path = env_path.as_posix()
    return Path(os.path.expandvars(env_path))


def rsync(src, tgt):
    os.system(f"rsync -rPvh {src} {tgt}")


def prsync(wild_src: str, tgt: str, ncores: int):
    os.system(f"ls -d {wild_src} | xargs -n 1 -P {ncores} -I% rsync -rvh % {tgt}")


def try_symlink(src: str, dst: str):
    """Create symlink if it doesn't exist

    Args:
        src (str): Source path
        dst (str): Destination path
    """
    # Create output dir if it doesn't exist
    try:
        os.symlink(src, dst)
        logger.info(f"Made symlink '{dst}'.")
    except FileExistsError:
        logger.info(f"Symlink '{dst}' exists.")


def try_mkdir(dir_path: str):
    """Create directory if it doesn't exist

    Args:
        dir_path (str): Path to directory
    """
    # Create output dir if it doesn't exist
    try:
        os.mkdir(dir_path)
        logger.info(f"Made directory '{dir_path}'.")
    except FileExistsError:
        logger.info(f"Directory '{dir_path}' exists.")


def gettable(tabledir: str, keyword: str) -> tuple[Table, str]:
    """Get a table from a directory given a keyword to glob.

    Args:
        tabledir (str): Directory.
        keyword (str): Keyword to glob for.
        verbose (bool, optional): Verbose output. Defaults to True.

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


def _samefile(src, dst):
    # Macintosh, Unix.
    if hasattr(os.path, "samefile"):
        try:
            return os.path.samefile(src, dst)
        except OSError:
            return False


def copyfile(src, dst, *, follow_symlinks=True, verbose=True):
    """Copy data from src to dst.

    If follow_symlinks is not set and src is a symbolic link, a new
    symlink will be created instead of copying the file it points to.

    """
    if _samefile(src, dst):
        raise SameFileError(f"{src!r} and {dst!r} are the same file")

    for fn in [src, dst]:
        try:
            st = os.stat(fn)
        except OSError:
            # File most likely does not exist
            pass
        else:
            # XXX What about other special files? (sockets, devices...)
            if stat.S_ISFIFO(st.st_mode):
                raise SpecialFileError(f"`{fn}` is a named pipe")

    if not follow_symlinks and os.path.islink(src):
        os.symlink(os.readlink(src), dst)
    else:
        with open(src, "rb") as fsrc:
            with open(dst, "wb") as fdst:
                copyfileobj(fsrc, fdst, verbose=verbose)
    return dst


def copyfileobj(fsrc, fdst, length=16 * 1024, verbose=True):
    # copied = 0
    total = os.fstat(fsrc.fileno()).st_size
    with tqdm(
        total=total,
        disable=(not verbose),
        unit_scale=True,
        desc="Copying file",
        file=TQDM_OUT,
    ) as pbar:
        while True:
            buf = fsrc.read(length)
            if not buf:
                break
            fdst.write(buf)
            copied = len(buf)
            pbar.update(copied)
