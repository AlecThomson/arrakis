#!/usr/bin/env python3

import os
from glob import glob
import tarfile
from tqdm.auto import tqdm
import dask
from dask import delayed
from dask.distributed import Client

@delayed
def tar_cubelets(out_dir: str, casda_dir: str, prefix: str) -> None:
    """Find and tar cubelets for a given source with a given prefix.

    Args:
        out_dir (str): Output directory
        casda_dir (str): CASDA directory containing cubelets/
        prefix (str): Prefix of cubelets to tar
    """
    print(f"Tarring {prefix}...")
    with tarfile.open(
        os.path.join(out_dir, f"{prefix}_cubelets.tar"), "w"
    ) as tar:
        _cube_list = glob(
            os.path.join(casda_dir, "cubelets", f"{prefix}*.fits")
        )
        for cube in _cube_list:
            tar.add(cube, arcname=os.path.basename(cube))
    print(f"...done {prefix}!")

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

    cube_list = glob(
        os.path.join(casda_dir, "cubelets", "*.fits")
    )
    print(f"{len(cube_list)} cublets to tar...")
    sources = set([os.path.basename(cube)[:13] for cube in tqdm(cube_list, desc="Sources")])
    print(f"...into {len(sources)} sources")
    out_dir = os.path.join(casda_dir, "cubelets_tar")
    os.makedirs(out_dir, exist_ok=True)
    print(f"Output directory: {out_dir}")
    outputs = []
    for source in tqdm(sources, desc="Tarring"):
        outputs.append(tar_cubelets(out_dir, casda_dir, source))

    dask.compute(*outputs)

    print("Done!")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("casda-dar", help="CASDA directory containing cublets/ to tar", type=str)
    args = parser.parse_args()
    # with Client() as client:
    main(args.cadsa_dir)
