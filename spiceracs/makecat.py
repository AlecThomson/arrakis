#!/usr/bin/env python3
"""Make a SPICE-RACS catalogue"""
import numpy as np
import warnings
from astropy.table import QTable, Column
from astropy.io import fits
import pymongo
from tqdm import tqdm, trange
from spiceracs import columns_possum
from spiceracs.utils import get_db, test_db
import rmtable.rmtable as RMT
import json
from IPython import embed


def main(
    field: str,
    host: str,
    username: str = None,
    password: str = None,
    verbose=True,
    outfile: str = None,
    cat_format: str = None,
) -> None:
    """Main

    Args:
        field (str): RACS field name
        host (str): MongoDB host IP
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo password. Defaults to None.
        verbose (bool, optional): Verbose output. Defaults to True.
        outfile (str, optional): Output file name. Defaults to None.
        cat_format (str, optional): Type of catalogue .e.g. fits. Defaults to None.
    """
    # default connection (ie, local)
    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )
    query = {
        "$and": [{f"beams.{field}": {"$exists": True}}, {f"beams.{field}.DR1": True}]
    }

    beams = beams_col.find(query).sort("Source_ID")
    all_island_ids = sorted(beams_col.distinct("Source_ID", query))
    query = {"$and": [{"Source_ID": {"$in": all_island_ids}}, {"rmclean1d": True}]}

    # tab = RMT.RMTable()
    tab = QTable()
    # Add items to main cat using RMtable standard
    for j, [name, typ, src, col, unit] in enumerate(
        tqdm(
            zip(
                columns_possum.output_cols,
                columns_possum.output_types,
                columns_possum.input_sources,
                columns_possum.input_names,
                columns_possum.output_units,
            ),
            total=len(columns_possum.output_cols),
            desc="Making table by column",
            disable=not verbose,
        ),
    ):
        data = []
        if src == "cat":
            for comp in comp_col.find(query):
                data += [comp[col]]
            new_col = Column(data=data, name=name, dtype=typ, unit=unit)
            tab.add_column(new_col)

        if src == "synth":
            for comp in comp_col.find(query):
                try:
                    data += [comp["rmsynth_summary"][col]]
                except KeyError:
                    data += [comp["rmclean_summary"][col]]
            new_col = Column(data=data, name=name, dtype=typ, unit=unit)
            tab.add_column(new_col)

        if src == "header":
            for comp in comp_col.find(query):
                data += [comp["header"][col]]
            new_col = Column(data=data, name=name, dtype=typ, unit=unit)
            tab.add_column(new_col)

    for selcol in tqdm(columns_possum.sourcefinder_columns, desc="Adding BDSF data"):
        data = []
        for comp in comp_col.find(query):
            data += [comp[selcol]]
        new_col = Column(data=data, name=selcol)
        tab.add_column(new_col)
    rmtab = RMT.from_table(tab)
    # Get Galatic coords
    rmtab["l"], rmtab["b"] = RMT.calculate_missing_coordinates_column(
        rmtab["ra"], rmtab["dec"], to_galactic=True
    )
    rmtab["rm_method"] = "RM Synthesis"
    rmtab["standard_telescope"] = "ASKAP"

    if outfile is None:
        print(rmtab)

    if outfile is not None:
        rmtab.table.write(outfile, format=cat_format, overwrite=True)

    if verbose:
        print("Done!")


def cli():
    """Command-line interface"""
    import argparse
    from astropy.utils.exceptions import AstropyWarning

    warnings.simplefilter("ignore", category=AstropyWarning)
    from astropy.io.fits.verify import VerifyWarning

    warnings.simplefilter("ignore", category=VerifyWarning)
    # Help string to be shown using the -h option
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
    SPICE-RACS Stage 7:
    Make RM catalogue.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field to mosaic - e.g. 2132-50A."
    )

    parser.add_argument(
        "host",
        metavar="host",
        type=str,
        help="Host of mongodb (probably $hostname -i).",
    )

    parser.add_argument(
        "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="verbose output [False]."
    )

    parser.add_argument(
        "-w",
        "--write",
        dest="outfile",
        default=None,
        type=str,
        help="File to save table to [None].",
    )

    parser.add_argument(
        "-f",
        "--format",
        dest="format",
        default=None,
        type=str,
        help="Format for output file [None].",
    )

    args = parser.parse_args()
    if args.outfile and not args.format:
        parser.error("Please provide an output file format.")

    verbose = args.verbose

    host = args.host
    test_db(
        host=args.host, username=args.username, password=args.password, verbose=verbose
    )

    main(
        field=args.field,
        host=host,
        username=args.username,
        password=args.password,
        verbose=verbose,
        limit=args.limit,
        outfile=args.outfile,
        cat_format=args.format,
    )


if __name__ == "__main__":
    cli()
