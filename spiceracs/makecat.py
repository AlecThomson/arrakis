import warnings
from astropy.table import QTable, Column
import pymongo
from tqdm import tqdm, trange
from spiceracs import columns_possum
import rmtable.rmtable as RMT


def main(args, verbose=False):
    """Main script.
    """
    outdir = args.outdir
    if outdir[-1] == '/':
        outdir = outdir[:-1]
    outdir = f'{outdir}/cutouts'
    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

    # Basic querey
    if args.pol and not args.unres:
        myquery = {"polarized": True}
        myquery = {"$and": [{"rmclean1d": True}, myquery]}
    elif args.unres and not args.pol:
        myquery = {"resolved": False}
        myquery = {"$and": [{"rmclean1d": True}, myquery]}
    elif args.pol and args.unres:
        myquery = {"$and": [{"rmclean1d": True}, {
            "resolved": False}, {"polarized": True}]}

    elif args.pol and not args.loners:
        myquery = {"polarized": True}
        myquery = {"$and": [{"rmclean1d": True}, myquery]}
    elif args.loners and not args.pol:
        myquery = {"n_components": 1}
        myquery = {"$and": [{"rmclean1d": True}, myquery]}
    elif args.pol and args.loners:
        myquery = {"$and": [{"rmclean1d": True}, {
            "n_components": 1}, {"polarized": True}]}

    else:
        myquery = {"rmclean1d": True}

    mydoc = mycol.find(myquery).sort("flux_peak", -1)
    count = mycol.count_documents(myquery)

    if args.limit is not None:
        count = args.limit

    #tab = RMT.RMTable()
    tab = QTable()

    # Add items to main cat using RMtable standard
    for j, [name, typ, src, col, unit] in tqdm(
            enumerate(
                zip(
                    columns_possum.output_cols,
                    columns_possum.output_types,
                    columns_possum.input_sources,
                    columns_possum.input_names,
                    columns_possum.output_units
                )
            ),
            desc='Making table',
            disable=not verbose
    ):
        data = []
        for i in range(count):
            try:
                for comp in range(mydoc[i]['n_components']):
                    if src == 'cat':
                        data.append(mydoc[i][f'component_{comp+1}'][col])
                    if src == 'synth': 
                        data.append(mydoc[i][f'comp_{comp+1}_rm_summary'][col]) 
                    if src == 'header': 
                        data.append(mydoc[i][f'header'][col]) 
                    else:
                        continue
            except KeyError:
                continue
        if data == []:
            continue
        new_col = Column(data=data, name=name, dtype=typ, unit=unit)
        tab.add_column(new_col)
    rmtab = RMT.from_table(tab)
    # Get Galatic coords
    rmtab['l'], rmtab['b'] = RMT.calculate_missing_coordinates_column(
        rmtab['ra'], rmtab['dec'],
        to_galactic=True
    )
    rmtab['rm_method'] = 'RM Synthesis'
    rmtab['standard_telescope'] = 'ASKAP'

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
        'outdir',
        metavar='outdir',
        type=str,
        help='Directory containing cutouts (in subdir outdir/cutouts).')

    parser.add_argument("-v", dest="verbose", action="store_true",
                        help="verbose output [False].")

    parser.add_argument("--pol", dest="pol", action="store_true",
                        help="Run on polarized sources [False].")

    parser.add_argument("--unres", dest="unres", action="store_true",
                        help="Run on unresolved sources [False].")

    parser.add_argument("--loners", dest="loners", action="store_true",
                        help="Run on single component sources [False].")

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

    if verbose:
        print('Testing MongoDB connection...')
    client = pymongo.MongoClient()  # default connection (ie, local)
    try:
        client.list_database_names()
    except pymongo.errors.ServerSelectionTimeoutError:
        raise Exception("Please ensure 'mongod' is running")
    else:
        if verbose:
            print('MongoDB connection succesful!')
    client.close()

    main(args, verbose=verbose)


if __name__ == "__main__":
    cli()
