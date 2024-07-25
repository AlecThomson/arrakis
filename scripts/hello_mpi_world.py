#!/usr/bin/env python
"""Parallel Hello World"""

from __future__ import annotations

import sys

from mpi4py import MPI


def main():
    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    name = MPI.Get_processor_name()

    sys.stdout.write(
        "Hello, World! I am process %d of %d on %s.\n" % (rank, size, name)
    )


def cli():
    import argparse

    parser = argparse.ArgumentParser(description="Run a parallel hello world")
    _ = parser.parse_args()
    main()


if __name__ == "__main__":
    cli()
