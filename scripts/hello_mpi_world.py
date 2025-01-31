#!/usr/bin/env python
"""
Parallel Hello World
"""

from __future__ import annotations

import sys

from mpi4py import MPI


def main():
    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    name = MPI.Get_processor_name()

    msg = f"Hello, World! I am process {rank} of {size} on {name}.\n"
    sys.stdout.write(msg)


def cli():
    import argparse

    parser = argparse.ArgumentParser(description="Run a parallel hello world")
    _ = parser.parse_args()
    main()


if __name__ == "__main__":
    cli()
