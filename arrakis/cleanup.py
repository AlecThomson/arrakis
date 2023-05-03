#!/usr/bin/env python3
"""DANGER ZONE: Purge directories of un-needed FITS files."""
import logging
import os
import time
from glob import glob
from typing import List, Union

from dask import delayed
from dask.distributed import Client, LocalCluster

from arrakis.logger import logger
from arrakis.utils import chunk_dask


@delayed
def cleanup(workdir: str, stoke: str) -> None:
    """Clean up beam images

    Args:
        workdir (str): Directory containing images
        stoke (str): Stokes parameter
    """
    # Clean up beam images
    # old_files = glob(f"{workdir}/*.cutout.*.{stoke.lower()}.*beam[00-36]*.fits")
    # for old in old_files:
    #     os.remove(old)

    pass


def main(
    datadir: str,
    stokeslist: Union[List[str], None] = None,
    verbose=True,
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

    if datadir is not None:
        if datadir[-1] == "/":
            datadir = datadir[:-1]

    cutdir = f"{datadir}/cutouts"
    files = sorted(
        [
            name
            for name in glob(f"{cutdir}/*")
            if os.path.isdir(os.path.join(cutdir, name))
        ]
    )

    outputs = []
    for file in files:
        if os.path.basename(file) == "slurmFiles":
            continue
        for stoke in stokeslist:
            output = cleanup(file, stoke)
            outputs.append(output)

    futures = chunk_dask(
        outputs=outputs,
        task_name="cleanup",
        progress_text="Running cleanup",
        verbose=verbose,
    )

    logger.info("Cleanup done!")


def cli():
    """Command-line interface"""
    import argparse

    logostr = """
     mmm   mmm   mmm   mmm   mmm
     )-(   )-(   )-(   )-(   )-(
    ( S ) ( P ) ( I ) ( C ) ( E )
    |   | |   | |   | |   | |   |
    |___| |___| |___| |___| |___|
     mmm     mmm     mmm     mmm
     )-(     )-(     )-(     )-(
    ( R )   ( A )   ( C )   ( S )
    |   |   |   |   |   |   |   |
    |___|   |___|   |___|   |___|

    """

    # Help string to be shown using the -h option
    descStr = f"""
    {logostr}
    Arrakis Stage:

    Clean up after LINMOS

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "outdir",
        metavar="outdir",
        type=str,
        help="Directory containing cutouts (in subdir outdir/cutouts).",
    )

    parser.add_argument(
        "-v", dest="verbose", action="store_true", help="Verbose output [False]."
    )
    args = parser.parse_args()
    verbose = args.verbose

    if verbose:
        logger.setLevel(logging.INFO)

    cluster = LocalCluster(n_workers=20)
    client = Client(cluster)

    main(datadir=args.outdir, stokeslist=None, verbose=verbose)

    client.close()
    cluster.close()


if __name__ == "__main__":
    cli()
