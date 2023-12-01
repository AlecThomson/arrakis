#!/usr/bin/env python3

import os
import tarfile
from glob import glob

import dask
from dask import delayed
from tqdm.auto import tqdm

from arrakis.logger import logger


@delayed
def tar_cubelets(out_dir: str, casda_dir: str, prefix: str) -> None:
    """Find and tar cubelets for a given source with a given prefix.

    Args:
        out_dir (str): Output directory
        casda_dir (str): CASDA directory containing cubelets/
        prefix (str): Prefix of cubelets to tar
    """
    logger.info(f"Tarring {prefix}...")
    with tarfile.open(os.path.join(out_dir, f"{prefix}_cubelets.tar"), "w") as tar:
        _cube_list = glob(os.path.join(casda_dir, "cubelets", f"{prefix}*.fits"))
        for cube in _cube_list:
            tar.add(cube, arcname=os.path.basename(cube))
    logger.info(f"...done {prefix}!")


def main(casda_dir: str):
    """Find cublets with unique prefixes and tar them.

    Args:
        casda_dir (str): CASDA directory containing cubelets/

    Raises:
        FileNotFoundError: If casda_dir does not exist or does not contain cubelets/
    """
    casda_dir = os.path.abspath(casda_dir)
    if not os.path.exists(casda_dir):
        raise FileNotFoundError(f"Directory {casda_dir} does not exist")
    if not os.path.exists(os.path.join(casda_dir, "cubelets")):
        raise FileNotFoundError(f"Directory {casda_dir} does not contain cubelets/")

    cube_list = glob(os.path.join(casda_dir, "cubelets", "*.fits"))
    logger.info(f"{len(cube_list)} cublets to tar...")
    sources = set(
        [os.path.basename(cube)[:13] for cube in tqdm(cube_list, desc="Sources")]
    )
    logger.info(f"...into {len(sources)} sources")
    out_dir = os.path.join(casda_dir, "cubelets_tar")
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f"Output directory: {out_dir}")
    outputs = []
    for source in tqdm(sources, desc="Tarring"):
        outputs.append(tar_cubelets(out_dir, casda_dir, source))

    dask.compute(*outputs)

    logger.info("Done!")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "casda_dir", help="CASDA directory containing cublets/ to tar", type=str
    )
    parser.add_argument(
        "-v", "--verbose", help="Increase output verbosity", action="store_true"
    )
    args = parser.parse_args()

    if args.verbose:
        logger.setLevel("INFO")

    main(args.casda_dir)
