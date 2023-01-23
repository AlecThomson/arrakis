#!/usr/bin/env python3

import os
from glob import glob
import tarfile
from tqdm.auto import tqdm
import dask
from dask import delayed
from dask.distributed import Client

@delayed
def tar_cubelets(out_dir, casda_dir, source):
    print(f"Tarring {source}...")
    with tarfile.open(
        os.path.join(out_dir, f"{source}_cubelets.tar"), "w"
    ) as tar:
        _cube_list = glob(
            os.path.join(casda_dir, "cubelets", f"{source}*.fits")
        )
        for cube in _cube_list:
            tar.add(cube, arcname=os.path.basename(cube))
    print(f"...done {source}!")

def main(thing):
    casda_dir = os.path.abspath(f"casda_{thing}")
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
    parser.add_argument("thing", help="Thing to tar")
    args = parser.parse_args()
    # with Client() as client:
    main(args.thing)
