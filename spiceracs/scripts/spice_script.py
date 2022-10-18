#!/usr/bin/env python3
"""Run spice scripts as subcommands"""
import argparse
import importlib
import os
import sys
from glob import glob


def main(argv=None):
    parser = argparse.ArgumentParser(prog="spice_script", add_help=True)
    here = os.path.dirname(os.path.abspath(__file__))
    globber = os.path.join(here, "*.py")
    scripts = glob(globber)
    script_names = [os.path.basename(s).split("/")[-1].split(".")[0] for s in scripts]
    badkeys = ("__init__", "spice_script")
    for key in badkeys:
        try:
            idx = script_names.index(key)
            script_names.pop(idx)
        except ValueError:
            pass
    parser.add_argument(choices=script_names, dest="subprogram", help="Subprogram to run")
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