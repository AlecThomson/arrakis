#!/usr/bin/env python3
"""DANGER ZONE: Purge directories of un-needed FITS files."""
import logging
import os
from glob import glob
from pathlib import Path
from typing import List, Union

from dask.distributed import Client, LocalCluster
from prefect import flow, task, unmapped

from arrakis.logger import logger
from arrakis.utils.pipeline import logo_str

logger.setLevel(logging.INFO)


@task(name="Cleanup directory")
def cleanup(workdir: str, stokeslist: List[str]) -> None:
    """Clean up beam images

    Args:
        workdir (str): Directory containing images
        stoke (str): Stokes parameter
    """
    if os.path.basename(workdir) == "slurmFiles":
        return
    for stoke in stokeslist:
        # Clean up beam images
        # old_files = glob(f"{workdir}/*.cutout.*.{stoke.lower()}.*beam[00-36]*.fits")
        # for old in old_files:
        #     os.remove(old)

        ...


@flow(name="Cleanup")
def main(
    datadir: Path,
    stokeslist: Union[List[str], None] = None,
) -> None:
    """Clean up beam images

    Args:
        datadir (str): Directory with sub dir 'cutouts'
        client (Client): Dask Client
        stokeslist (List[str], optional): List of Stokes parameters to purge. Defaults to None.
        verbose (bool, optional): Verbose output. Defaults to True.
    """
    if stokeslist is None:
        stokeslist = ["I", "Q", "U", "V"]

    cutdir = datadir / "cutouts"
    files = sorted(
        [
            name
            for name in glob(f"{cutdir}/*")
            if os.path.isdir(os.path.join(cutdir, name))
        ]
    )
    outputs = cleanup.map(
        workdir=files,
        stokeslist=unmapped(stokeslist),
    )

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

    main(datadir=Path(args.outdir), stokeslist=None, verbose=verbose)


if __name__ == "__main__":
    cli()
