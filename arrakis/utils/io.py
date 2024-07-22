"""I/O utilities"""

from __future__ import annotations

import logging
import os
import shlex
import subprocess as sp
import warnings
from pathlib import Path

from astropy.utils.exceptions import AstropyWarning
from spectral_cube.utils import SpectralCubeWarning

from arrakis.logger import TqdmToLogger, logger
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
        verbose (bool, optional): Verbose output. Defaults to True.
    """
    # Create output dir if it doesn't exist
    try:
        os.symlink(src, dst)
        logger.info(f"Made symlink '{dst}'.")
    except FileExistsError:
        logger.info(f"Symlink '{dst}' exists.")
