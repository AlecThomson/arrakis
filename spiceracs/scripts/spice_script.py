#!/usr/bin/env python3
"""Run spice scripts as subcommands"""
import argparse
import importlib
import os
import sys
from glob import glob


def main(argv=None):
    parser = argparse.ArgumentParser(prog="spice_script", add_help=True)
    scripts = glob(os.path.join(__file__,"*.py"))
    scripts = [s.split("/")[-1].split(".")[0] for s in scripts]
    # Pop the __init__.py
    try:
        idx = scripts.index("__init__")
        scripts.pop(idx)
    except ValueError:
        pass
    parser.add_argument(choices=scripts, dest="subprogram", help="Subprogram to run")
    opts = parser.parse_args(argv)
    return opts.subprogram

def cli():
    # Get the subprogram
    subprogram = main([sys.argv[1]])
    module = importlib.import_module(f"spiceracs.scripts.{subprogram}")
    module.cli()

# Initialize
if __name__ == "__main__":
    cli()