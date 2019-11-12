import numpy as np
from glob import glob
from spectral_cube import SpectralCube
from tqdm import trange, tqdm
from astropy.wcs import WCS
from astropy.table import Table
import sys
import os


def getdata(cubedir, tabledir, verbose=True):
    """Get the spectral and source-finding data.

    Args:
        cubedir: Directory containing data cubes in FITS format.
        tabledir: Directory containing Selavy results.

    Kwargs:
        verbose (bool): Whether to print messages.

    Returns:
        datadict (dict): Dictionary of necessary astropy tables and
            Spectral cubes.

    """
    # Glob out the necessary files
    # Data cubes
    cubes = glob(f'{cubedir}/image.restored.*contcube*linmos*fits')
    selavyfits = glob(f'{tabledir}/comp*.fits')  # Selavy images
    voisle = glob(f'{tabledir}/*island*.xml')  # Selvay VOTab

    if verbose:
        print('Getting islands from:', voisle, '\n')
        print('Getting spectral data from:', cubes, '\n')
        print('Getting source location data from:', selavyfits[0], '\n')

    # Get selvay data from VOTab
    i_tab = Table.read(voisle[0], format='votable')

    # Read data using Spectral cube
    i_taylor = SpectralCube.read(selavyfits[0], mode='denywrite')
    wcs_taylor = WCS(i_taylor.header)
    i_cube = SpectralCube.read(cubes[0], mode='denywrite')
    wcs_cube = WCS(i_cube.header)
    q_cube = SpectralCube.read(cubes[1], mode='denywrite')
    u_cube = SpectralCube.read(cubes[2], mode='denywrite')

    datadict = {
        "i_tab": i_tab,
        "i_taylor": i_taylor,
        "wcs_taylor": wcs_taylor,
        "wcs_cube": wcs_cube,
        "i_cube": i_cube,
        "q_cube": q_cube,
        "u_cube": u_cube
    }

    return datadict


def makecutout(datadict, outdir='.', pad=0, dryrun=False, verbose=True):

    # Get bounding boxes in WCS
    ra_min, dec_min, freq = datadict['wcs_taylor'].all_pix2world(
        datadict['i_tab']['col_x_min'], datadict['i_tab']['col_y_min'], 0, 0)
    ra_max, dec_max, freq = datadict['wcs_taylor'].all_pix2world(
        datadict['i_tab']['col_x_max'], datadict['i_tab']['col_y_max'], 0, 0)

    # Get bounding boxes in cube pixels
    x_min, y_min, _ = np.array(datadict['wcs_cube'].all_world2pix(
        ra_min, dec_min, freq, 0)).astype(int)
    x_max, y_max, _ = np.array(datadict['wcs_cube'].all_world2pix(
        ra_max, dec_max, freq, 0)).astype(int)
    dy, dx = y_max - y_min, x_max-x_min

    # Init cutouts
    i_cutouts = []
    q_cutouts = []
    u_cutouts = []
    outdir = f'{outdir}/cutouts'
    if dryrun:
        print('Dry run -- not saving to disk.')
    else:
        print(f'Saving to {outdir}/')
    source_dict_list = []

    i_cube = datadict['i_cube']
    q_cube = datadict['q_cube']
    u_cube = datadict['u_cube']
    if not dryrun:
        try:
            os.mkdir(outdir)
            print('Made directory.')
        except FileExistsError:
            print('Directory exists.')

    # TO-DO: Max cut on size
    for i in trange(
        len(datadict['i_tab']),
        disable=(not verbose),
        desc='Extracting cubelets'
    ):
        source_dict = {}
        # Skip if source is outside of cube bounds
        if (y_max[i] > i_cube.shape[2] or x_max[i] > i_cube.shape[2] or
                x_min[i] < 0 or y_min[i] < 0):
            continue
        # Check if pad puts bbox outside of cube
        elif (int(y_min[i]-pad*dy[i]) > 0 and
              int(x_min[i]-pad*dx[i]) > 0 and
              int(y_max[i]+pad*dy[i]) < i_cube.shape[2] and
              int(x_max[i]+pad*dx[i]) < i_cube.shape[2]):

            i_cutout = i_cube[
                :,
                int(y_min[i]-pad*dy[i]):int(y_max[i]+pad*dy[i]),
                int(x_min[i]-pad*dx[i]):int(x_max[i]+pad*dx[i])
            ]
            q_cutout = q_cube[
                :,
                int(y_min[i]-pad*dy[i]):int(y_max[i]+pad*dy[i]),
                int(x_min[i]-pad*dx[i]):int(x_max[i]+pad*dx[i])
            ]
            u_cutout = u_cube[
                :,
                int(y_min[i]-pad*dy[i]):int(y_max[i]+pad*dy[i]),
                int(x_min[i]-pad*dx[i]):int(x_max[i]+pad*dx[i])
            ]

            i_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
            q_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
            u_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
            i_cutouts.append(i_cutout)
            q_cutouts.append(q_cutout)
            u_cutouts.append(u_cutout)

            source_dict['header'] = i_cutout.header
            for name in datadict['i_tab'].colnames:
                source_dict[name.replace('col_', '')
                            ] = datadict['i_tab'][name][i]
            source_dict_list.append(source_dict)
        else:
            i_cutout = i_cube[:, y_min[i]:y_max[i], x_min[i]:x_max[i]]
            q_cutout = q_cube[:, y_min[i]:y_max[i], x_min[i]:x_max[i]]
            u_cutout = u_cube[:, y_min[i]:y_max[i], x_min[i]:x_max[i]]
            i_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
            q_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
            u_cutout.meta['OBJECT'] = datadict['i_tab']['col_island_name'][i]
            i_cutouts.append(i_cutout)
            q_cutouts.append(q_cutout)
            u_cutouts.append(u_cutout)

            source_dict['header'] = i_cutout.header
            for name in datadict['i_tab'].colnames:
                source_dict[name.replace('col_', '')
                            ] = datadict['i_tab'][name][i]
            source_dict_list.append(source_dict)

    # Set up locations where files will be saved
    for i in trange(
        len(source_dict_list),
        disable=(not verbose),
        desc='Finding locations'
    ):
        for stoke in ['i', 'q', 'u']:
            name = source_dict_list[i]['island_name']
            outname = f'{outdir}/{name}.cutout.{stoke}.fits'
            source_dict_list[i][f'{stoke}_file'] = outname

    cutouts = {
        "i": i_cutouts,
        "q": q_cutouts,
        "u": u_cutouts
    }
    return cutouts, source_dict_list


def getbytes(cutout):
    return cutout[0, :, :].nbytes*cutout.shape[0]*1e-6


def getsize(pool, cutouts):
    if args.mpi:
        sizes_bytes = list(
            pool.map(getbytes, [cutouts[i] for i in range(len(cutouts))])
        )
    else:
        sizes_bytes = list(tqdm(
            pool.imap_unordered(getbytes, [cutouts[i]
                                           for i in range(len(cutouts))]),
            total=len(cutouts),
            desc='Getting sizes')
        )
    sizes_bytes = np.array(sizes_bytes)
    print('Size in MB: ', sizes_bytes.sum())
    print('Size in GB: ', sizes_bytes.sum()/1000)


def writefits(arg):
    """Write cutouts to disk.

    Writes a cutout, as stored in source_dict, to disk. The file
    location should already be specified in source_dict. This format is
    intended for parallel use with pool.map syntax.

    Args:
        arg: The tuple of (i, stoke)
            i: The index of the source in source_dict to write to disk.
                stoke: Which Stokes to write. Is a string of either 'i',
                    'q', or 'u'.
    """
    source_dict, cutout, stoke = arg
    outfile = source_dict[f'{stoke}_file']
    cutout.write(outfile, format='fits', overwrite=True)


def writeloop(pool, cutouts, source_dict_list):
    for stoke in ['i', 'q', 'u']:
        if args.mpi:
            pool.map(
                writefits,
                [[source_dict_list[i], cutouts[stoke][i], stoke]
                    for i in range(len(cutouts[stoke]))]
            )
        else:
            list(tqdm(
                pool.imap_unordered(
                    writefits,
                    [[source_dict_list[i], cutouts[stoke][i], stoke]
                     for i in range(len(cutouts[stoke]))]
                ),
                total=len(cutouts[stoke]),
                desc=f'Stokes {stoke}''
            )
            )


def main(pool=None, args=None, verbose=True):
    """The main script.
    """
    # Sort out args
    cubedir = args.cubedir
    tabledir = args.tabledir
    if cubedir[-1] == '/':
        cubedir = cubedir[:-1]

    if tabledir[-1] == '/':
        tabledir = tabledir[:-1]

    # Read in data
    if verbose:
        print('Reading data...')
    datadict = getdata(cubedir, tabledir, verbose=verbose)

    # Make cutouts
    outdir = args.outdir
    if outdir[-1] == '/':
        outdir = outdir[:-1]
    pad = args.pad
    dryrun = args.dryrun
    if verbose:
        print('Making cutouts....')
    cutouts, source_dict_list = makecutout(
        datadict,
        outdir=outdir,
        pad=pad,
        dryrun=dryrun,
        verbose=verbose
    )

    if args.getsize:
        if verbose:
            print('Checking size of single cube...')
        getsize(pool, cutouts['i'])

    if not dryrun:
        if verbose:
            print('Writing to disk...')
        writeloop(pool, cutouts, source_dict_list)

    pool.close()
    if verbose:
        print('Done!')


if __name__ == "__main__":
    import argparse
    import schwimmbad

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
        'cubedir',
        metavar='cubedir',
        type=str,
        help='Directory containing data cubes in FITS format.')

    parser.add_argument(
        'tabledir',
        metavar='tabledir',
        type=str,
        help='Directory containing Selavy results.')

    parser.add_argument(
        'outdir',
        metavar='outdir',
        type=str,
        help='Directory to store cutouts.')

    parser.add_argument(
        'pad',
        metavar='pad',
        type=float,
        default=0,
        help='Fractional padding around islands [0 -- no padding].')

    parser.add_argument(
        "-v",
        dest="verbose",
        action="store_true",
        help="Verbose output [False]."
    )

    parser.add_argument(
        "-d",
        dest="dryrun",
        action="store_true",
        help="Do a dry-run [False]."
    )

    parser.add_argument(
        "-s",
        dest="getsize",
        action="store_true",
        help="Estimate size of cutouts [False]."
    )

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        "--ncores",
        dest="n_cores",
        default=1,
        type=int, help="Number of processes (uses multiprocessing)."
    )
    group.add_argument(
        "--mpi",
        dest="mpi",
        default=False,
        action="store_true",
        help="Run with MPI."
    )

    args = parser.parse_args()
    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)
    verbose = args.verbose

    if args.mpi:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)

    main(pool=pool, args=args, verbose=verbose)
