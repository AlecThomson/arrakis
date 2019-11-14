#!/usr/bin/env python
import numpy as np
from astropy.table import Table
import pymongo
from tqdm import tqdm, trange
import warnings
from glob import glob
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pdb


def getdata(tabledir, verbose=True):
    """Get the spectral and source-finding data.

    Args:
        tabledir: Directory containing Selavy results.

    Kwargs:
        verbose (bool): Whether to print messages.

    Returns:
        datadict (dict): Dictionary of necessary astropy tables and
            Spectral cubes.

    """
    # Glob out the necessary files
    vocomp = glob(f'{tabledir}/*components*.xml')  # Selvay VOTab

    if verbose:
        print('Getting components from:', vocomp[0], '\n')

    # Get selvay data from VOTab
    components = Table.read(vocomp[0], format='votable')

    return components


def updatedb_comp(mycol, loners, verbose=True):
    """Update MongoDB with component data.

    The Selavy Island catalogue have unreliable fluxes, so we need to
    obtain the component fluxes. 

    Args:
        mycol (collection): MongoDB collection to match and update.
        loners (Table): Single-component sources with their Island_ID
            indexed.

    Kwargs:
        verbose (bool): Whether to print messages.

    """
    # Table columns to grab
    fields = ['col_flux_peak', 'col_flux_peak_err',
              'col_flux_int', 'col_flux_int_err']

    # Init loop
    skipped = 0
    updated = 0
    for i in trange(
        len(loners),
        disable=(not verbose),
        desc='Adding component values'
    ):
        # Find single-component sources that match Island ID
        myquery = {'$and': [{'n_components': {'$lt': 2}},
                            {"island_id":  loners['col_island_id'][i]}]}
        count = mycol.count_documents(myquery)

        if count == 0:
            skipped += 1
        else:
            updated += 1
        for field in fields:
            newvalues = {"$set": {field.replace('col', 'comp'):
                                  float(loners[field][i])}}
            result = mycol.update_one(myquery, newvalues)

    if verbose:
        print(f'{updated} sources updated.')
        print(f'{skipped} sources skipped.')


def getflux(mycol, verbose=True):
    """Get flux values.

    Args:
        mycol (collection): MongoDB collection to get data from.

    Kwargs:
        verbose (bool): Whether to print messages.

    Returns:
        flux_dict (dict): Arrays of peak flux, itegrated flux, and noise
            from source components.

    """
    # Find single component sources
    myquery = {'$and': [{"n_components": {"$lt": 2}},
                        {'comp_flux_peak': {"$exists": True}}]}
    mydoc = mycol.find(myquery) 
    count = mycol.count_documents(myquery)
    if verbose:
        print(f'Using {count} sources.')
    f_peak = []
    f_int = []
    noise = []
    for x in tqdm(
        mydoc,
        disable=(not verbose),
        desc='Getting flux and noise',
        total=count

    ):
        f_peak.append(x['comp_flux_peak'])
        noise.append(x['background_noise'])
        f_int.append(x['comp_flux_int'])

    f_peak = np.array(f_peak)
    noise = np.array(noise)
    f_int = np.array(f_int)
    snr = f_peak/noise
    flux_dict = {
        'peak': f_peak,
        'noise': noise,
        'int': f_int,
        'snr': snr
    }
    return flux_dict


def makeplot(flux_dict, snrcut, verbose=True):
    """Produce flux plot.

    Args:
        flux_dict (dict): Arrays of peak flux, itegrated flux, and noise
            from source components.

    Kwargs:
        verbose (bool): Whether to print messages.

    """
    warnings.filterwarnings("ignore",
                            message="divide by zero encountered in true_divide")

    plt.figure()
    plt.scatter(flux_dict['peak'], flux_dict['peak']/(flux_dict['int']),
                c=flux_dict['snr'], marker='.', norm=LogNorm(),
                vmin=snrcut, vmax=np.nanmax(flux_dict['snr'])
                )
    cbar = plt.colorbar()
    cbar.set_label('SNR ($S_{peak}$/background noise)')
    plt.xscale('log')
    plt.xlabel(r'S$_{peak}$ [mJy/beam]')
    plt.ylabel(r'S$_{peak}$/S$_{int}$')

    x = np.linspace(1e0, flux_dict['peak'].max(), num=100000)
    y = 1 - (0.1/np.log10(x))
    plt.plot(x, y, label='$1 - 0.1/log(S_{peak})$')

    x = np.linspace(1e0, flux_dict['peak'].max(), num=100000)
    y = 1 + (0.1/np.log10(x))

    plt.plot(x, y, label='$1 + 0.1/log(S_{peak})$')
    plt.ylim(-0.1, 1.5)
    plt.xlim(7e-1, 2e3)

    plt.legend()
    # Save plot
    outfig = 'fluxplot.png'
    plt.savefig(outfig, dpi=1000)
    if verbose:
        print(f'Plot saved to {outfig}.')


def updatedb_ur(mycol, snrcut, verbose=True):
    """Update MongoDB with unresolved tag.

    Args:
        mycol (collection): MongoDB collection to match and update.

    Kwargs:
        verbose (bool): Whether to print messages.

    """
    myquery = {'$and': [{"n_components": {"$lt": 2}},
                        {'comp_flux_peak': {"$exists": True}}]}
    mydoc = mycol.find(myquery)
    count = mycol.count_documents(myquery)
    if verbose:
        print(f'Updating {count} sources...')

    for item in tqdm(
        mydoc,
        disable=(not verbose),
        desc='Updating resolved sources',
        total=count
    ):
        # Get component fluxes
        s_peak = item['comp_flux_peak']
        s_int = item['comp_flux_int']
        noise = item['background_noise']
        snr = s_peak/noise
        m_id = item['_id']

        # Banfield (2015) cut
        test = 1 - (0.1 / np.log10(s_peak))
        frac = s_peak / s_int
        query = {"_id": m_id}
        if frac < test and snr > snrcut:
            newvalues = {"$set": {'resolved':  True}}
            mycol.update_one(query, newvalues)
        elif snr > snrcut:
            newvalues = {"$set": {'resolved':  False}}
            mycol.update_one(query, newvalues)
        # except:
        #    continue
    myquery = {"resolved": True}
    count = mycol.count_documents(myquery)
    if verbose:
        print(f'Identified {count} resolved sources.')
    myquery = {"resolved": False}
    count = mycol.count_documents(myquery)
    if verbose:
        print(f'Identified {count} unresolved sources.')


def main(args, verbose=True):
    """Main script.
    """
    tabledir = args.tabledir
    if tabledir[-1] == '/':
        tabledir = tabledir[:-1]

    if verbose:
        print('Reading catalogue...')
    components = getdata(tabledir, verbose=verbose)

    # Find single-component sources
    if verbose:
        print('Finding single-component sources...')
    loners = components[components['col_has_siblings'] == 0]
    loners.add_index('col_island_id')
    if verbose:
        print(f'Found {len(loners)} single-component sources.')

    if not args.database:
        if verbose:
            print('DB update required to continue.')
        return

    else:
        # Init DB
        client = pymongo.MongoClient()  # default connection (ie, local)
        mydb = client['racs']  # Create/open database
        mycol = mydb['spice']  # Create/open collection

        # Update DB
        if verbose:
            print('Adding component data to DB...')
        updatedb_comp(mycol, loners, verbose=verbose)

        # Retrieve data
        if verbose:
            print('Getting data...')
        flux_dict = getflux(mycol, verbose=verbose)

        if args.do_plot:
            if verbose:
                print('Making plot...')
            makeplot(flux_dict, args.snrcut, verbose=verbose)

        # Update DB
        if verbose:
            print('Adding resolved data to DB...')
            print(f'Using an SNR cut of {args.snrcut}.')
        updatedb_ur(mycol, args.snrcut, verbose=verbose)

        client.close()
        if verbose:
            print('Done!')


if __name__ == "__main__":
    import argparse
    import schwimmbad
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)

    # Help string to be shown using the -h option
    descStr = """
    SPICE-RACS Stage 2:
    Find unresolved sources
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        'tabledir',
        metavar='tabledir',
        type=str,
        help='Directory containing Selavy results.')

    parser.add_argument(
        "--snr",
        dest="snrcut",
        type=float,
        default=10,
        help="Signal-to-noise cut [10].")

    parser.add_argument(
        "-v",
        dest="verbose",
        action="store_true",
        help="Verbose output [False]."
    )

    parser.add_argument(
        "-m",
        dest="database",
        action="store_true",
        help="Add data to MongoDB [False]."
    )

    parser.add_argument(
        "-p",
        dest="do_plot",
        action="store_true",
        help="Show plots [False]."
    )

    args = parser.parse_args()
    verbose = args.verbose

    if args.database:
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
