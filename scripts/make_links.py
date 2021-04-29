#!/usr/bin/env python
import os
import subprocess
from shlex import split
from glob import glob
import argparse

def main(indir, outdir):
    images = glob(f"{os.path.abspath(indir)}/image.restored.*.contcube.*.fits")
    weights = glob(f"{os.path.abspath(indir)}/weights.*.contcube.*.fits")
    for f in images:
        name = os.path.basename(f)
        link = name.replace('.fits','.conv.fits')
        cmd = f"ln -s {f} {os.path.abspath(outdir)}/{link}"
        print(cmd)
        subprocess.run(split(cmd))

    for f in weights:
        name = os.path.basename(f)
        cmd = f"ln -s {f} {os.path.abspath(outdir)}/{name}"
        print(cmd)
        subprocess.run(split(cmd))

if __name__ == "__main__":
    descStr = f"""
    Create symlinks to ASKAP cubes in one directory to another.
    """
    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
    'indir',
    type=str,
    help='Directory containing cubes'
    )


    parser.add_argument(
    'outdir',
    type=str,
    help='Target output directory'
    )
    
    args = parser.parse_args()

    main(args.indir, args.outdir)