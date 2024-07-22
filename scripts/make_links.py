#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import subprocess
from glob import glob
from shlex import split

from arrakis.logger import logger, logging

logger.setLevel(logging.INFO)


def main(indir, outdir):
    images = glob(f"{os.path.abspath(indir)}/image.restored.*.contcube.*.fits")
    weights = glob(f"{os.path.abspath(indir)}/weights.*.contcube.*.fits")
    for f in images:
        name = os.path.basename(f)
        link = name.replace(".fits", ".conv.fits")
        cmd = f"ln -s {f} {os.path.abspath(outdir)}/{link}"
        logger.info(cmd)
        subprocess.run(split(cmd), check=False)

    for f in weights:
        name = os.path.basename(f)
        cmd = f"ln -s {f} {os.path.abspath(outdir)}/{name}"
        logger.info(cmd)
        subprocess.run(split(cmd), check=False)


def cli():
    descStr = """
    Create symlinks to ASKAP cubes in one directory to another.
    """
    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("indir", type=str, help="Directory containing cubes")

    parser.add_argument("outdir", type=str, help="Target output directory")

    args = parser.parse_args()

    main(args.indir, args.outdir)


if __name__ == "__main__":
    cli()
