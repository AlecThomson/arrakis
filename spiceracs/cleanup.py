#!/usr/bin/env python3
"""DANGER ZONE: Purge directories of un-needed FITS files."""
import os
from glob import glob
import time
from spiceracs.utils import chunk_dask
from dask import delayed
from dask.distributed import Client, LocalCluster
from typing import List
import logging as log


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
    datadir: str, client: Client, stokeslist: List[str] = None, verbose=True
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
        client=client,
        task_name="cleanup",
        progress_text="Running cleanup",
        verbose=verbose,
    )

    log.info("Cleanup done!")


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
    SPICE-RACS Stage:

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
        log.basicConfig(
            level=log.INFO,
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            force=True
        )
    else:
        log.basicConfig(
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            force=True
        )

    cluster = LocalCluster(n_workers=20)
    client = Client(cluster)

    main(datadir=args.outdir, client=client, stokeslist=None, verbose=verbose)

    client.close()
    cluster.close()


if __name__ == "__main__":
    cli()
