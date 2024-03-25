#!/usr/bin/env python3
"""DANGER ZONE: Purge directories of un-needed FITS files."""
import logging
import tarfile
from pathlib import Path
from typing import List, Union

from prefect import flow, task, unmapped
from tqdm.auto import tqdm

from arrakis.logger import TqdmToLogger, logger
from arrakis.utils.pipeline import logo_str

logger.setLevel(logging.INFO)

TQDM_OUT = TqdmToLogger(logger, level=logging.INFO)


@task(name="Purge cublets")
def purge_cubelet_beams(filepath: Path) -> Path:
    """Clean up beam images

    Args:
        workdir (str): Directory containing images
        stoke (str): Stokes parameter
    """

    # Clean up beam images
    logger.critical(f"Removing {filepath}")
    filepath.unlink()
    logger.info(f"Beam image purged from {filepath}")
    return filepath


@task(name="Make cutout tarball")
def make_cutout_tarball(cutdir: Path) -> Path:
    """Make a tarball of the cutouts directory

    Args:
        cutdir (Path): Directory containing cutouts

    Returns:
        Path: Path to the tarball
    """
    logger.info(f"Making tarball of {cutdir}")
    tarball = cutdir.with_suffix(".tar.gz")
    with tarfile.open(tarball, "w:gz") as tar:
        tar.add(cutdir, arcname=cutdir.name)
    return tarball


@flow(name="Cleanup")
def main(
    datadir: Path,
) -> None:
    """Clean up beam images

    Args:
        datadir (Path): Directory with sub dir 'cutouts'
    """

    cutdir = datadir / "cutouts"

    # First, make a tarball of the cutouts
    tarball = make_cutout_tarball.submit(cutdir=cutdir)
    logger.info(f"Tarball created: {tarball.result()}")

    # Purge the beam images from cutouts
    to_purge_cutouts = cutdir.glob("*/*.beam[00-36]*.fits")
    logger.info(f"Purging beam images from {cutdir}")
    purged = purge_cubelet_beams.map(to_purge_cutouts)
    logger.info(f"Beam images purged: {len([p.result() for p in purged])}")

    # Purge the big beams
    to_purge_bigbeams = datadir.glob("*.beam[00-36]*.fits")
    logger.info(f"Purging big beams from {datadir}")
    purged_bigbeams = purge_cubelet_beams.map(to_purge_bigbeams)
    logger.info(f"Big beams purged: {len([p.result() for p in purged_bigbeams])}")

    logger.info("Cleanup done!")


def cli():
    """Command-line interface"""
    import argparse

    # Help string to be shown using the -h option
    descStr = f"""
    {logo_str}
    Arrakis Stage:

    Clean up after LINMOS

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "outdir",
        metavar="outdir",
        type=Path,
        help="Directory containing cutouts (in subdir outdir/cutouts).",
    )

    parser.add_argument(
        "-v", dest="verbose", action="store_true", help="Verbose output [False]."
    )
    args = parser.parse_args()
    verbose = args.verbose

    if verbose:
        logger.setLevel(logging.DEBUG)

    main(datadir=Path(args.outdir))


if __name__ == "__main__":
    cli()
