#!/usr/bin/env python
import numpy as np
from astropy.table import Table
import pymongo
from tqdm import tqdm, trange
import warnings
from glob import glob
import matplotlib.pyplot as plt

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
    vocomp = glob(f'{tabledir}/*components*.xml') # Selvay VOTab

    if verbose:
        print('Getting components from:', vocomp[0], '\n')

    # Get selvay data from VOTab
    components = Table.read(vocomp[0], format='votable')

    return components

def updatedb(mycol, loners, verbose=True):
    fields = ['col_flux_peak', 'col_flux_peak_err', 'col_flux_int', 'col_flux_int_err']
    for i in trange(
            len(loners),
            disable=(not verbose),
            desc='Adding component values'
        ):
        myquery = { '$and': [ { 'n_components': { '$lt': 2 } }, \
                        { "island_id":  loners['col_island_id'][i] } ] }

        for field in fields:
            newvalues = { "$set": { field.replace('col','comp'):  float(loners[field][i])} }
            mycol.update_one(myquery, newvalues)

def getflux(mycol, verbose=True):
        """
        """
        myquery = { "n_components": { "$lt": 2 } }
        mydoc = mycol.find(myquery)
        count = mydoc.count_documents()
        if verbose: print(f'Using {count} sources.')
        f_peak = []
        f_int = []
        noise = []
        for x in tqdm(
                mydoc,
                disable=(not verbose),
                desc='Getting flux and noise'
            ):
            try:
                f_peak.append(x['comp_flux_peak'])
                noise.append(x['background_noise'])
                f_int.append(x['comp_flux_int'])
            except KeyError:
                continue
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

def makeplot(flux_dict):
    plt.figure()
    plt.scatter(f_peak, f_peak/(f_int), c=snr, marker='.', norm=LogNorm())
    cbar = plt.colorbar()
    cbar.set_label('SNR ($S_{peak}$/background noise)')
    plt.xscale('log')
    plt.xlabel(r'S$_{peak}$ [mJy/beam]')
    plt.ylabel(r'S$_{peak}$/S$_{int}$')

    x = np.linspace(1e0, f_peak.max(), num=100000)
    y = 1 - (0.1/np.log10(x))

    plt.plot(x, y, label='$1 - 0.1/log(S_{peak})$')


    x = np.linspace(1e0, f_peak.max(), num=100000)
    y = 1 + (0.1/np.log10(x))

    plt.plot(x, y, label='$1 + 0.1/log(S_{peak})$')
    plt.ylim(-0.1,1.5)
    plt.xlim(7e-1, 2e3)

    plt.legend()

def main(args, verbose=True):
    """Main script.
    """
    tabledir = args.tabledir
    if tabledir[-1] == '/':
        tabledir = tabledir[:-1]

    if verbose: print('Reading catalogue...')
    components = getdata(tabledir, verbose=verbose)


    # Find single-component sources
    if verbose: print('Finding single-component sources...')
    loners = components[components['col_has_siblings']==0]
    loners.add_index('col_island_id')
    if verbose: print(f'Found {len(loners)} single-component sources.')

    if not args.database:
        if verbose: print('DB update required to continue.')
        return
    
    else:
        # Init DB
        client = pymongo.MongoClient()  # default connection (ie, local)
        mydb = client['racs']  # Create/open database
        mycol = mydb['spice']  # Create/open collection

        # Update DB
        if verbose: print('Adding component data to DB...')
        updatedb(mycol, loners, verbose=verbose)

        # Retrieve data
        if verbose: print('Getting data...')
        getflux(flux_dict)


if __name__ == "__main__":
    import argparse
    import schwimmbad
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)

    # Help string to be shown using the -h option
    descStr = """
    Produce cutouts of a given RACS field.
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


    args = parser.parse_args()
    verbose = args.verbose

    if args.database:
        if verbose: print('Testing MongoDB connection...')
        client = pymongo.MongoClient()  # default connection (ie, local)
        try:
            client.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")
        else:
            if verbose: print('MongoDB connection succesful!')
        client.close()

    main(args, verbose=verbose)