import numpy as np
import warnings
from astropy.table import QTable, Column
from astropy.io import fits
import pymongo
from tqdm import tqdm, trange
from spiceracs import columns_possum
import rmtable.rmtable as RMT
import json
from IPython import embed


def main(args, verbose=False):
    """Main script.
    """
    outdir = args.outdir
    if outdir[-1] == '/':
        outdir = outdir[:-1]
    outdir = f'{outdir}/cutouts'
    field = args.field
    host = args.host
    # default connection (ie, local)
    with pymongo.MongoClient(host=host) as client:
        mydb = client['spiceracs']  # Create/open database
        isl_col = mydb['islands']  # Create/open collection
        comp_col = mydb['components']  # Create/open collection
        beams_col = mydb['beams']  # Create/open collection

    query = {
        '$and':  [
            {f'beams.{field}': {'$exists': True}},
            {f'beams.{field}.DR1': True}
        ]
    }

    beams = beams_col.find(query).sort('Source_ID')
    all_island_ids = sorted(beams_col.distinct('Source_ID', query))
    query = {
        '$and': [
            {
                'Source_ID': {'$in': all_island_ids}
            },
            {
                'rmclean1d': True
            }
        ]
    }
    count = comp_col.count_documents(query)

    if args.limit is not None:
        count = args.limit

    #tab = RMT.RMTable()
    tab = QTable()
    # Add items to main cat using RMtable standard
    for j, [name, typ, src, col, unit] in enumerate(
        tqdm(
            zip(
                columns_possum.output_cols,
                columns_possum.output_types,
                columns_possum.input_sources,
                columns_possum.input_names,
                columns_possum.output_units
            ),
            total=len(columns_possum.output_cols),
            desc='Making table by column',
            disable=not verbose),
    ):
        data = []
        if src == 'cat':
            for comp in comp_col.find(query):
                data += [comp[col]]
            new_col = Column(data=data, name=name, dtype=typ, unit=unit)
            tab.add_column(new_col)

        if src == 'synth':
            for comp in comp_col.find(query):
                try:
                    data += [comp['rmsynth_summary'][col]]
                except KeyError:
                        data += [comp['rmclean_summary'][col]]
            new_col = Column(data=data, name=name, dtype=typ, unit=unit)
            tab.add_column(new_col)

        if src == 'header':
            for comp in comp_col.find(query):
                data += [comp['header'][col]]
            new_col = Column(data=data, name=name, dtype=typ, unit=unit)
            tab.add_column(new_col)

    for selcol in tqdm(columns_possum.sourcefinder_columns, desc='Adding BDSF data'):
        data = []
        for comp in comp_col.find(query):
            data += [comp[selcol]]
        new_col = Column(data=data, name=selcol)
        tab.add_column(new_col)
    rmtab = RMT.from_table(tab)
    # Get Galatic coords
    rmtab['l'], rmtab['b'] = RMT.calculate_missing_coordinates_column(
        rmtab['ra'], rmtab['dec'],
        to_galactic=True
    )
    rmtab['rm_method'] = 'RM Synthesis'
    rmtab['standard_telescope'] = 'ASKAP'

    if args.outfile is None:
        print(rmtab)

    if args.outfile is not None:
        rmtab.table.write(args.outfile, format=args.format, overwrite=True)

    if verbose:
        print('Done!')

def cli():
    """Command-line interface
    """
    import argparse
    import schwimmbad
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
    from astropy.io.fits.verify import VerifyWarning
    warnings.simplefilter('ignore', category=VerifyWarning)
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
    parser = argparse.ArgumentParser(description=descStr,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "field",
        metavar="field",
        type=str,
        help="RACS field to mosaic - e.g. 2132-50A."
    )
    parser.add_argument(
        'outdir',
        metavar='outdir',
        type=str,
        help='Directory containing cutouts (in subdir outdir/cutouts).')

    parser.add_argument(
        'host',
        metavar='host',
        type=str,
        help='Host of mongodb (probably $hostname -i).')

    parser.add_argument("-v", dest="verbose", action="store_true",
                        help="verbose output [False].")

    parser.add_argument("--limit", dest="limit", default=None,
                        type=int, help="Limit number of sources [All].")

    parser.add_argument("-w", "--write", dest="outfile", default=None,
                        type=str, help="File to save table to [None].")

    parser.add_argument("-f", "--format", dest="format", default=None,
                        type=str, help="Format for output file [None].")

    args = parser.parse_args()
    if args.outfile and not args.format:
        parser.error('Please provide an output file format.')

    verbose = args.verbose

    host = args.host
    if verbose:
        print('Testing MongoDB connection...')
    # default connection (ie, local)
    with pymongo.MongoClient(host=host) as client:
        try:
            client.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")
        else:
            if verbose:
                print('MongoDB connection succesful!')

    main(args, verbose=verbose)


if __name__ == "__main__":
    cli()
