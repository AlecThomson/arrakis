#!/usr/bin/env python
import os
import stat
from tqdm import tqdm
from glob import glob
from spectral_cube import SpectralCube
from astropy.table import Table
from astropy.wcs import WCS
import astropy.units as u
import functools
print = functools.partial(print, flush=True)


def getfreq(cube, outdir=None, filename=None, verbose=True):
    """Get list of frequencies from FITS data.

    Gets the frequency list from a given cube. Can optionally save
    frequency list to disk.

    Args:
        cube (str or SpectralCube): File or cube to get spectral
            axis from. If a file, it will be opened using SpectralCube.

    Kwargs:
        outdir (str): Where to save the output file. If not given, data
            will not be saved to disk.

        filename (str): Name of frequency list file. Requires 'outdir'
            to also be specified.

        verbose (bool): Whether to print messages.

    Returns:
        freq (list): Frequencies of each channel in the input cube.

    """

    # If cube is a file, open with SpectralCube
    if type(cube) is str:
        cube = SpectralCube.read(cube, mode='denywrite')

    # Test that cube is Spectral cube
    assert cube is SpectralCube, "cube should be a SpectralCube!"

    # Get frequencies
    freq = cube.spectral_axis

    # Write to file if outdir is specified
    if outdir is not None:
        if outdir[-1] == '/':
            outdir = outdir[:-1]
        if filename is None:
            outfile = f'{outdir}/frequencies.txt'
        else:
            outfile = f'{outdir}/{filename}'
        if verbose:
            print(f'Saving to {outfile}')
        np.savetext(outfile, freq)
    return freq


def gettable(tabledir, keyword, verbose=True):
    """Get the spectral and source-finding data.

    Args:
        tabledir (str): Directory containing Selavy results.
        keyword (str): Glob out files containing '*.keyword.*'.

    Kwargs:
        verbose (bool): Whether to print messages.

    Returns:
        datadict (dict): Dictionary of necessary astropy tables and
            Spectral cubes.

    """
    if tabledir[-1] == '/':
        tabledir = tabledir[:-1]
    # Glob out the necessary files
    files = glob(f'{tabledir}/*.{keyword}*.xml')  # Selvay VOTab
    filename = files[0]
    if verbose:
        print(f'Getting table data from {filename}...')

    # Get selvay data from VOTab
    table = Table.read(filename, format='votable')

    return table, filename


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
    if cubedir[-1] == '/':
        cubedir = cubedir[:-1]

    if tabledir[-1] == '/':
        tabledir = tabledir[:-1]
    # Glob out the necessary files
    # Data cubes
    cubes = glob(f'{cubedir}/image.restored.*contcube*linmos.fits')
    selavyfits = glob(f'{tabledir}/comp*.fits')  # Selavy images
    # Get selvay data from VOTab
    i_tab, voisle = gettable(
        tabledir, 'islands', verbose=verbose)  # Selvay VOTab

    if verbose:
        print(f'Getting spectral data from: {cubes}', '\n')
        print(f'Getting source location data from:', selavyfits[0], '\n')

    # Read data using Spectral cube
    i_taylor = SpectralCube.read(selavyfits[0], mode='denywrite')
    wcs_taylor = WCS(i_taylor.header)
    i_cube = SpectralCube.read(cubes[0], mode='denywrite')
    wcs_cube = WCS(i_cube.header)
    q_cube = SpectralCube.read(cubes[1], mode='denywrite')
    u_cube = SpectralCube.read(cubes[2], mode='denywrite')

    # Mask out using Stokes I == 0 -- seems to be the current fill value
    mask = ~(i_cube == 0*u.jansky/u.beam)
    i_cube = i_cube.with_mask(mask)
    q_cube = q_cube.with_mask(mask)
    u_cube = u_cube.with_mask(mask)

    datadict = {
        "i_tab": i_tab,
        "i_taylor": i_taylor,
        "wcs_taylor": wcs_taylor,
        "wcs_cube": wcs_cube,
        "i_cube": i_cube,
        "q_cube": q_cube,
        "u_cube": u_cube,
        "i_file": cubes[0],
        "q_file": cubes[1],
        "u_file": cubes[2],
    }

    return datadict


class Error(OSError):
    pass


class SameFileError(Error):
    """Raised when source and destination are the same file."""


class SpecialFileError(OSError):
    """Raised when trying to do a kind of operation (e.g. copying) which is
    not supported on a special file (e.g. a named pipe)"""


class ExecError(OSError):
    """Raised when a command could not be executed"""


class ReadError(OSError):
    """Raised when an archive cannot be read"""


class RegistryError(Exception):
    """Raised when a registry operation with the archiving
    and unpacking registeries fails"""


def _samefile(src, dst):
    # Macintosh, Unix.
    if hasattr(os.path, 'samefile'):
        try:
            return os.path.samefile(src, dst)
        except OSError:
            return False


def copyfile(src, dst, *, follow_symlinks=True, verbose=True):
    """Copy data from src to dst.

    If follow_symlinks is not set and src is a symbolic link, a new
    symlink will be created instead of copying the file it points to.

    """
    if _samefile(src, dst):
        raise SameFileError("{!r} and {!r} are the same file".format(src, dst))

    for fn in [src, dst]:
        try:
            st = os.stat(fn)
        except OSError:
            # File most likely does not exist
            pass
        else:
            # XXX What about other special files? (sockets, devices...)
            if stat.S_ISFIFO(st.st_mode):
                raise SpecialFileError("`%s` is a named pipe" % fn)

    if not follow_symlinks and os.path.islink(src):
        os.symlink(os.readlink(src), dst)
    else:
        with open(src, 'rb') as fsrc:
            with open(dst, 'wb') as fdst:
                copyfileobj(fsrc, fdst, verbose=verbose)
    return dst


def copyfileobj(fsrc, fdst, length=16*1024, verbose=True):
    #copied = 0
    total = os.fstat(fsrc.fileno()).st_size
    with tqdm(
            total=total,
            disable=(not verbose),
            unit_scale=True,
            desc='Copying file'
    ) as pbar:
        while True:
            buf = fsrc.read(length)
            if not buf:
                break
            fdst.write(buf)
            copied = len(buf)
            pbar.update(copied)
