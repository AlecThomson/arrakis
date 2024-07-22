#!/usr/bin/env python3
"""DANGER ZONE: Purge directories of un-needed FITS files."""

import argparse
import logging
import shutil
import tarfile
from pathlib import Path
from typing import List

import astropy.units as u
import numpy as np
from prefect import flow, get_run_logger, task
from tqdm.auto import tqdm

from arrakis.logger import TqdmToLogger, UltimateHelpFormatter, logger
from arrakis.utils.io import verify_tarball
from arrakis.utils.pipeline import generic_parser, logo_str

logger.setLevel(logging.INFO)

TQDM_OUT = TqdmToLogger(logger, level=logging.INFO)


# @task(name="Purge cublets")
def purge_cubelet_beams(filepath: Path) -> Path:
    """Clean up beam images

    Args:
        workdir (str): Directory containing images
        stoke (str): Stokes parameter
    """
    # Clean up beam images
    logger.critical(f"Removing {filepath}")
    filepath.unlink(missing_ok=True)
    logger.info(f"Beam image purged from {filepath}")
    return filepath


@task(name="Make cutout tarball")
def make_cutout_tarball(cutdir: Path, overwrite: bool = False) -> Path:
    """Make a tarball of the cutouts directory

    Args:
        cutdir (Path): Directory containing cutouts

    Returns:
        Path: Path to the tarball
    """
    logger = get_run_logger()
    logger.info(f"Making tarball of {cutdir}")
    tarball = cutdir.with_suffix(".tar")
    if tarball.exists() and not overwrite:
        logger.warning(f"Tarball {tarball} exists. Skipping.")
        return tarball
    all_things = list(cutdir.glob("*"))
    with tarfile.open(tarball, "w") as tar:
        for cutout in tqdm(all_things, file=TQDM_OUT, desc="Tarballing cutouts"):
            tar.add(cutout, arcname=cutout.name)

    logger.info(f"Tarball created: {tarball}")

    verification = verify_tarball(tarball)
    if not verification:
        raise RuntimeError(f"Verification of {tarball} failed!")

    logger.critical(f"Removing {cutdir}")
    shutil.rmtree(cutdir)
    return tarball


@flow(name="Cleanup")
def main(
    datadir: Path,
    overwrite: bool = False,
) -> None:
    """Clean up beam images flow

    Args:
        datadir (Path): Directory with sub dir 'cutouts'
        overwrite (bool): Overwrite existing tarball
    """

    cutdir = datadir / "cutouts"

    # First, make a tarball of the cutouts
    tarball = make_cutout_tarball.submit(cutdir=cutdir, overwrite=overwrite)
    logger.info(f"Tarball created: {tarball.result()}")

    # Purge the big beams
    to_purge_bigbeams = list(datadir.glob("*.beam[00-36]*.fits"))
    logger.warning(f"Will purge {len(to_purge_bigbeams)} big beams")
    to_purge_pickles = list(datadir.glob("*.pkl"))
    logger.warning(f"Will purge {len(to_purge_pickles)} beam pickles")
    to_purge_weights = list(datadir.glob("weights*.txt"))
    logger.warning(f"Will purge {len(to_purge_weights)} weight files")

    # Purge the cublet beams
    to_purge_cublets = list(cutdir.glob("*/*.beam[00-36]*.fits"))
    logger.warning(f"Will purge {len(to_purge_cublets)} cublet beams")

    to_purge_all = (
        to_purge_bigbeams + to_purge_cublets + to_purge_pickles + to_purge_weights
    )

    total_file_size = np.sum([p.stat().st_size for p in to_purge_all]) * u.byte
    logger.warning(f"Purging {len(to_purge_all)} files from {datadir}")
    logger.warning(f"Will free {total_file_size.to(u.GB)}")
    purged: List[Path] = []
    for to_purge in tqdm(to_purge_all, file=TQDM_OUT, desc="Purging big beams"):
        purged.append(purge_cubelet_beams(to_purge))
    logger.info(f"Files purged: {len(purged)}")

    logger.info("Cleanup done!")


def cleanup_parser(parent_parser: bool = False) -> argparse.ArgumentParser:
    # Help string to be shown using the -h option
    descStr = f"""
    {logo_str}
    Arrakis Stage:

    Clean up after LINMOS

    """

    # Parse the command line options
    cleanup_parser = argparse.ArgumentParser(
        add_help=not parent_parser,
        description=descStr,
        formatter_class=UltimateHelpFormatter,
    )
    parser = cleanup_parser.add_argument_group("cleanup arguments")
    parser.add_argument(
        "--overwrite",
        dest="overwrite",
        action="store_true",
        help="Overwrite existing tarball",
    )

    return cleanup_parser


def cli():
    """Command-line interface"""
    gen_parser = generic_parser(parent_parser=True)
    clean_parser = cleanup_parser(parent_parser=True)
    parser = argparse.ArgumentParser(
        parents=[gen_parser, clean_parser],
        formatter_class=UltimateHelpFormatter,
        description=clean_parser.description,
    )
    args = parser.parse_args()

    verbose = args.verbose

    if verbose:
        logger.setLevel(logging.DEBUG)

    main(datadir=Path(args.datadir), overwrite=args.overwrite)


if __name__ == "__main__":
    cli()
