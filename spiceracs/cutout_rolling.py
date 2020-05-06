#!/usr/bin/env python

from astropy.table import Table, vstack
from glob import glob
import os
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm, trange


def source2beams(ra, dec, database, max_sep=4):
    """Find RACS beams containing a position.

    Arguments:
        ra {float} -- RA of source in degrees.
        dec {[type]} -- DEC of source in degrees.
        database {astropy.table.table.Table} -- RACS database loaded as one Table.
    Keyword Arguments:
        max_sep {int} -- Maximum angular distance to beam centre in degrees (default: {4})

    Returns:
        beams {astropy.table.table.Table} -- RACS database rows matching the source location.
    """
    c1 = SkyCoord(database['RA_DEG']*u.deg,
                  database['DEC_DEG']*u.deg, frame='icrs')
    c2 = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    sep = c1.separation(c2)
    beams = database[sep < max_sep*u.deg]
    return beams


def main(args):
    """Main script

    Arguments:
        args {[type]} -- commandline args
    """
    # Get data base
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    basedir = f"{scriptdir}/../askap_surveys/RACS/admin/epoch_0"
    beamfiles = glob(f"{basedir}/beam_inf*")

    # Init first field
    beamfile = beamfiles[0]
    database = Table.read(beamfile)
    basename = os.path.basename(beamfile)
    idx = basename.find('RACS_test')
    FIELD = basename[idx:-4]
    SBID = basename[9:idx-1]
    database.add_column(FIELD, name='FIELD_NAME', index=0)
    database.add_column(int(SBID), name='SBID', index=0)

    # Add in all others
    for i, beamfile in enumerate(tqdm(beamfiles, desc='Reading RACS database')):
        if i == 0:
            continue
        else:
            tab = Table.read(beamfile)
            basename = os.path.basename(beamfile)
            idx = basename.find('RACS_test')
            FIELD = basename[idx:-4]
            SBID = basename[9:idx-1]
            tab.add_column(FIELD, name='FIELD_NAME', index=0)
            tab.add_column(int(SBID), name='SBID', index=0)
            database = vstack([database, tab])
    
    # Read in main catalogue
    island_cat = Table.read(
        '/group/askap/chale/RACS/RACS-25asec/Selavy-Mosaiced/Final-Catalogues/RACS-25asec-Mosaiced_Islands_NoArtefacts_TF_v2020_04_22.fits')

idx = island_cat['Pointing_ID'] == 'RACS_test4_1.05_2132-50A'
island_cat = island_cat[idx]
beam_dict = {}
for row in tqdm(island_cat):
    beam_dict.update(
        {row['island_id']: source2beams(row['ra_deg_cont'], row['dec_deg_cont'], database, 1)}
    )
    
    
def cli():
    """Command-line interface
    """
    import argparse

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

    descStr = f"""
    {logostr}
    SPICE-RACS Stage 1:
    Produce cubelets from a RACS field using a Selavy table.
    If Stokes V is present, it will be squished into RMS spectra.

    To use with MPI:
       mpirun -n $NPROCS python -u cutout.py $cubedir $tabledir
       $outdir --mpi
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )

    args = parser.parse_args()

    main(args)


if __name__ == "__main__":
    cli()
