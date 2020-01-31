import rmtable.rmtable as RMT
import pymongo
import sys
import warnings
from tqdm import tqdm, trange


def main(args, verbose=False):
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
    if verbose:
        print('Done!')

    tab = RMT.RMTable()

    # Add items to main cat using RMtable standard
    for i in trange(
        count,
        desc='Generating RMtable',
        total=count,
        disable=(not verbose)
    ):
        tab.table.add_row(
            {
                'ra': mydoc[i]['ra_deg_cont'],
                'dec': mydoc[i]['dec_deg_cont'],
                'rm': mydoc[i]['rm_summary']['phiPeakPIchan_rm2'],
                'rm_err': mydoc[i]['rm_summary']['dPhiPeakPIchan_rm2']
            }
        )

    # Get Galatic coords
    tab['l'], tab['b'] = RMT.calculate_missing_coordinates_column(
        tab['ra'], tab['dec'],
        to_galactic=True
    )

    tab.write_tsv('/Users/tho822/Desktop/junk.fits')

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


    args = parser.parse_args()

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
