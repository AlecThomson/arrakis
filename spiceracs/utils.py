#!/usr/bin/env python
import numpy as np
import os
import stat
from tqdm import tqdm
from glob import glob
from spectral_cube import SpectralCube
import subprocess
from pathlib import Path
from astropy.table import Table
from astropy.wcs import WCS
import astropy.units as u
import functools
print = functools.partial(print, flush=True)


def cpu_to_use(max_cpu, count):
    """Find number of cpus to use.

    Find the right number of cpus to use when dividing up a task, such
    that there are no remainders.

    Args:
        max_cpu (int): Maximum number of cores to use for a process.
        count (float): Number of tasks.
    
    Returns:
        Maximum number of cores to be used that divides into the number
        of tasks (int).
    """
    factors = []
    for i in range(1, count + 1):
        if count % i == 0:
            factors.append(i)
    factors = np.array(factors)
    return max(factors[factors <= max_cpu])


def tmatchtwo(inN, valuesN, matcher='sky', params=10, omode='out',
              out='tmatch.default.xml', join='1or2', verbose=True):
    """
    inN = <tableN>       (StarTable)
        The location of input table #N. This may take one of the
        following forms:
            A filename.
            A URL.
            The special value "-", meaning standard input. In this case
            the input format must be given explicitly using the ifmtN
            parameter. Note that not all formats can be streamed in this
            way.
            A system command line with either a "<" character at the
            start, or a "|" character at the end ("<syscmd" or
            "syscmd|"). This executes the given pipeline and reads from
            its standard output. This will probably only work on
            unix-like systems.

    valuesN = <expr-list>       (String[])
        Defines the values from table N which are used to determine
        whether a match has occurred. These will typically be coordinate
        values such as RA and Dec and perhaps some per-row error values
        as well, though exactly what values are required is determined
        by the kind of match as determined by matcher. Depending on the
        kind of match, the number and type of the values required will
        be different. Multiple values should be separated by whitespace;
        if whitespace occurs within a single value it must be 'quoted'
        or "quoted". Elements of the expression list are commonly just
        column names, but may be algebraic expressions calculated from
        zero or more columns as explained in Section 10.

    matcher = <matcher-name>       (MatchEngine)
        Defines the nature of the matching that will be performed. 
        Depending on the name supplied, this may be positional matching
        using celestial or Cartesian coordinates, exact matching on the
        value of a string column, or other things. A list and
        explanation of the available matching algorithms is given in
        Section 7.1. The value supplied for this parameter determines
        the meanings of the values required by the params, values* and
        tuning parameter(s).
        [Default: sky]

    params = <match-params>       (String[])
        Determines the parameters of this match. This is typically one
        or more tolerances such as error radii. It may contain zero or
        more values; the values that are required depend on the match
        type selected by the matcher parameter. If it contains multiple 
        alues, they must be separated by spaces; values which contain a
        space can be 'quoted' or "quoted".    

    omode = out|meta|stats|count|cgi|discard|topcat|samp|plastic
                |tosql|gui       (ProcessingMode)
        The mode in which the result table will be output. The default
        mode is out, which means that the result will be written as a
        new table to disk or elsewhere, as determined by the out and
        ofmt parameters. However, there are other possibilities, which
        correspond to uses to which a table can be put other than
        outputting it, such as displaying metadata, calculating
        statistics, or populating a table in an SQL database. For some
        values of this parameter, additional parameters (<mode-args>)
        are required to determine the exact behaviour.
        [Default: out]

    out = <out-table>       (TableConsumer)
        The location of the output table. This is usually a filename to
        write to. If it is equal to the special value "-" (the default)
        the output table will be written to standard output.
        This parameter must only be given if omode has its default
        value of "out".
        [Default: -]

    join = 1and2|1or2|all1|all2|1not2|2not1|1xor2       (JoinType)
        Determines which rows are included in the output table. The
        matching algorithm determines which of the rows from the first
        table correspond to which rows from the second. This parameter
        determines what to do with that information. Perhaps the most
        obvious thing is to write out a table containing only rows which
        correspond to a row in both of the two input tables. However,
        you may also want to see the unmatched rows from one or both
        input tables, or rows present in one table but unmatched in the
        other, or other possibilities. The options are:
            1and2: An output row for each row represented in both input
                tables (INNER JOIN)
            1or2: An output row for each row represented in either or
                both of the input tables (FULL OUTER JOIN)
            all1: An output row for each matched or unmatched row in
                table 1 (LEFT OUTER JOIN)
            all2: An output row for each matched or unmatched row in
                table 2 (RIGHT OUTER JOIN)
            1not2: An output row only for rows which appear in the first
                table but are not matched in the second table
            2not1: An output row only for rows which appear in the
                second table but are not matched in the first table
            1xor2: An output row only for rows represented in one of the
                input tables but not the other one
        [Default: 1and2]

    """
    assert len(inN) == 2, 'Can only match 2 tables!'
    assert len(valuesN) == 2, 'Can only match 2 tables!'
    assert len(inN) == len(valuesN), 'Need same number of inputs (2).'
    stiltspath = Path(os.path.realpath(__file__)
                      ).parent.parent/"thirdparty"/"stilts"/"stilts.jar"
    # Construct command
    if verbose:
        progress = 'log'
    if not verbose:
        progress = 'none'
    command = ['java', '-jar', stiltspath, 'tmatch2',
               f'matcher={matcher}', f'params={params}',
               f'omode={omode}', f'out={out}', f'progress={progress}',
               f'join={join}'
               ]
    for i in range(len(inN)):
        command.append(f'in{i+1}={inN[i]}')
        command.append(f'values{i+1}={valuesN[i]}')
    command = [str(i) for i in command]
    # Run STILTS
    proc = subprocess.run(command,
                          capture_output=(not verbose),
                          encoding="utf-8", check=True)


def tmatchn(nin, inN, valuesN,
            matcher='sky', params=10, omode='out',
            out='tmatch.default.xml', verbose=True):
    """Run STILTS tmatchn
    nin = <count>       (Integer)
        The number of input tables for this task. For each of the input
        tables N there will be associated parameters ifmtN, inN and
        icmdN.

    inN = <tableN>       (StarTable)
        The location of input table #N. This may take one of the
        following forms:
            A filename.
            A URL.
            The special value "-", meaning standard input. In this case
            the input format must be given explicitly using the ifmtN
            parameter. Note that not all formats can be streamed in this
            way.
            A system command line with either a "<" character at the
            start, or a "|" character at the end ("<syscmd" or
            "syscmd|"). This executes the given pipeline and reads from
            its standard output. This will probably only work on
            unix-like systems.

    valuesN = <expr-list>       (String[])
        Defines the values from table N which are used to determine
        whether a match has occurred. These will typically be coordinate
        values such as RA and Dec and perhaps some per-row error values
        as well, though exactly what values are required is determined
        by the kind of match as determined by matcher. Depending on the
        kind of match, the number and type of the values required will
        be different. Multiple values should be separated by whitespace;
        if whitespace occurs within a single value it must be 'quoted'
        or "quoted". Elements of the expression list are commonly just
        column names, but may be algebraic expressions calculated from
        zero or more columns as explained in Section 10.

    matcher = <matcher-name>       (MatchEngine)
        Defines the nature of the matching that will be performed.
        Depending on the name supplied, this may be positional matching
        using celestial or Cartesian coordinates, exact matching on the
        value of a string column, or other things. A list and
        explanation of the available matching algorithms is given in
        Section 7.1. The value supplied for this parameter determines
        the meanings of the values required by the params, values* and
        tuning parameter(s).
        [Default: sky]

    params = <match-params>       (String[])
        Determines the parameters of this match. This is typically one
        or more tolerances such as error radii. It may contain zero or
        more values; the values that are required depend on the match
        type selected by the matcher parameter. If it contains multiple
        values, they must be separated by spaces; values which contain a
        space can be 'quoted' or "quoted".

    omode = out|meta|stats|count|cgi|discard|topcat|samp|plastic|tosql|gui
            (ProcessingMode)
        The mode in which the result table will be output. The default
        mode is out, which means that the result will be written as a
        new table to disk or elsewhere, as determined by the out and
        ofmt parameters. However, there are other possibilities, which
        correspond to uses to which a table can be put other than
        outputting it, such as displaying metadata, calculating
        statistics, or populating a table in an SQL database. For some
        values of this parameter, additional parameters (<mode-args>)
        are required to determine the exact behaviour.
        [Default: out]
    out = <out-table>       (TableConsumer)
        The location of the output table. This is usually a filename to
        write to. If it is equal to the special value "-" (the default)
        the output table will be written to standard output.
        This parameter must only be given if omode has its default value 
        of "out".
        [Default: -]
    """
    stiltspath = Path(os.path.realpath(__file__)
                      ).parent.parent/"thirdparty"/"stilts"/"stilts.jar"
    # Construct command
    if verbose:
        progress = 'log'
    if not verbose:
        progress = 'none'
    command = ['java', '-jar', stiltspath, 'tmatchn', f'nin={nin}',
               f'matcher={matcher}', f'params={params}',
               f'omode={omode}', f'out={out}', f'progress={progress}']
    for i in range(len(inN)):
        command.append(f'in{i+1}={inN[i]}')
        command.append(f'values{i+1}={valuesN[i]}')
    command = [str(i) for i in command]
    # Run STILTS
    proc = subprocess.run(command,
                          capture_output=(not verbose),
                          encoding="utf-8", check=True)


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
    assert type(cube) is SpectralCube, "cube should be a SpectralCube!"

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


def getdata(cubedir='./', tabledir='./', mapdir='./', verbose=True):
    """Get the spectral and source-finding data.

    Args:
        cubedir: Directory containing data cubes in FITS format.
        tabledir: Directory containing Selavy results.
        mapdir: Directory containing MFS image.

    Kwargs:
        verbose (bool): Whether to print messages.

    Returns:
        datadict (dict): Dictionary of necessary astropy tables and
            Spectral cubes.

    """
    if cubedir[-1] == '/':
        cubedir = cubedir[:-1]

    if mapdir[-1] == '/':
        mapdir = mapdir[:-1]

    if tabledir[-1] == '/':
        tabledir = tabledir[:-1]
    # Glob out the necessary files
    # Data cubes
    icubes = glob(f'{cubedir}/image.restored.i.*contcube*linmos.fits')
    qcubes = glob(f'{cubedir}/image.restored.q.*contcube*linmos.fits')
    ucubes = glob(f'{cubedir}/image.restored.u.*contcube*linmos.fits')
    vcubes = glob(f'{cubedir}/image.restored.v.*contcube*linmos.fits')

    cubes = [icubes, qcubes, ucubes, vcubes]
    # Selavy images
    selavyfits = glob(f'{mapdir}/image.i*.linmos.taylor.0.restored.fits')
    # Get selvay data from VOTab
    i_tab, voisle = gettable(
        tabledir, 'islands', verbose=verbose)  # Selvay VOTab

    if verbose:
        print(f'Getting spectral data from: {cubes}', '\n')
        print(f'Getting source location data from:', selavyfits[0], '\n')

    # Read data using Spectral cube
    i_taylor = SpectralCube.read(selavyfits[0], mode='denywrite')
    wcs_taylor = WCS(i_taylor.header)
    i_cube = SpectralCube.read(icubes[0], mode='denywrite')
    wcs_cube = WCS(i_cube.header)
    q_cube = SpectralCube.read(qcubes[0], mode='denywrite')
    u_cube = SpectralCube.read(ucubes[0], mode='denywrite')

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
        "i_file": icubes[0],
        "q_file": qcubes[0],
        "u_file": ucubes[0],
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
