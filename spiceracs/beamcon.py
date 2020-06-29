#!/usr/bin/env python
from spectral_cube.utils import SpectralCubeWarning
import warnings
import os
import stat
import sys
import numpy as np
import scipy.signal
from astropy import units as u
from astropy.io import fits
from spectral_cube import SpectralCube
from radio_beam import Beam, Beams
from radio_beam.utils import BeamError
import schwimmbad
from tqdm import tqdm, trange
from IPython import embed
from glob import glob
import au2
import functools
from functools import partial
import psutil
#print = functools.partial(print, f'[{psutil.Process().cpu_num()}]', flush=True)
from mpi4py import MPI
mpiComm = MPI.COMM_WORLD
n_cores = mpiComm.Get_size()
print = functools.partial(print, f'[{mpiComm.rank}]', flush=True)
warnings.filterwarnings(action='ignore', category=SpectralCubeWarning,
                        append=True)

#############################################
#### ADAPTED FROM SCRIPT BY T. VERNSTROM ####
#############################################


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


def round_up(n, decimals=0):
    multiplier = 10 ** decimals
    return np.ceil(n * multiplier) / multiplier


def getbeams(beamdir, verbose=False):
    """

    colnames=['Channel', 'BMAJarcsec', 'BMINarcsec', 'BPAdeg']
    """

    stokes = ["i", "q", "u", "v"]
    beams = [str(n).zfill(2) for n in range(36)]
    beamdict = {}
    for stoke in stokes:
        beamdict.update(
            {
                stoke: {}
            }
        )
        for beam in beams:
            beamlog = glob(f"{beamdir}/beamlog*.{stoke}.*beam{beam}.txt")[0]
            beamdata = np.genfromtxt(beamlog, names=True)
            nchan = beamdata.shape[0]
            beamdict[stoke].update(
                {
                    beam: {
                        'beamlog': beamlog,
                        'beams': beamdata,
                        'nchan': nchan
                    }
                }
            )
    return beamdict


def getfacs(oldbeams, conbeams, masks, dx, dy, verbose=False):
    """Get beam info
    """
    facs = []
    for oldbeam, conbeam, mask in zip(oldbeams, conbeams, masks):
        if mask:
            facs.append(np.nan)
        else:
            if verbose:
                print(f"Current beam is", oldbeam)
            fac, amp, outbmaj, outbmin, outbpa = au2.gauss_factor(
                [
                    conbeam.major.to(u.arcsec).value,
                    conbeam.minor.to(u.arcsec).value,
                    conbeam.pa.to(u.deg).value
                ],
                beamOrig=[
                    oldbeam.major.to(u.arcsec).value,
                    oldbeam.minor.to(u.arcsec).value,
                    oldbeam.pa.to(u.deg).value
                ],
                dx1=dx.to(u.arcsec).value,
                dy1=dy.to(u.arcsec).value
            )
            facs.append(fac)
    return facs


def smooth(image, grid, conbeam, sfactor, verbose=False):
    """Do the smoothing
    """
    if np.isnan(sfactor):
        return image*np.nan
    if np.isnan(image).all():
        return image
    else:
        # using Beams package
        if verbose:
            print(f'Using convolving beam', conbeam)
        gauss_kern = conbeam.as_kernel(grid)

        conbm1 = gauss_kern.array/gauss_kern.array.max()
        newim = scipy.signal.convolve(
            image.astype('f8'), conbm1, mode='same')
    newim *= sfactor
    return newim


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


def worker(file, beamdict, target_beam, dryrun=True, verbose=False):
    with fits.open(file, memmap=True, mode='denywrite') as hdulist:
        header = hdulist[0].header
        fitscube = np.squeeze(hdulist[0].data)
    dxas = header['CDELT1']*-1*u.deg
    dyas = header['CDELT2']*u.deg
    # Get beam info
    dirname = os.path.dirname(file)
    filename = os.path.basename(file)
    basename = filename[filename.find('image.'):]
    if dirname == '':
        dirname = '.'
    stoke = basename[15+basename.find('image.restored.')                     :15+basename.find('image.restored.')+1]
    beamno = basename[basename.find(
        'beam')+4:basename.find('beam')+4+2]
    beamlog = beamdict[stoke][beamno]['beamlog']

    beams = beamdict[stoke][beamno]['PSFs']
    nchan = beamdict[stoke][beamno]['nchan']
    conbeams = beamdict[stoke][beamno]['conbeams']
    mask = beamdict[stoke][beamno]['mask']

    facs = getfacs(beams, conbeams, ~mask, dxas, dyas, verbose=False)

    newcube = np.ones_like(fitscube) * np.nan
    if verbose:
        print(f'Smoothing {filename}')
    for i, (plane, conbeam, sfactor) in enumerate(
        zip(fitscube, conbeams, facs)
    ):
        newim = smooth(plane, dyas, conbeam, sfactor, verbose=False)
        newcube[i] = newim

    if not dryrun:
        outname = "sm." + filename
        outfile = f'{dirname}/{outname}'
        if verbose:
            print(f"Saved to {outfile}")
        header = target_beam.attach_to_header(header)
        fits.writeto(outfile, newcube, header, overwrite=True)


def main(pool, args, verbose=True):
    if args.dryrun:
        if verbose:
            print('Doing a dry run -- no files will be saved')
    # Fix up outdir
    logdir = args.logdir
    if logdir is not None:
        if logdir[-1] == '/':
            logdir = logdir[:-1]

    beamdict = getbeams(logdir)

    # Get data for all beams
    bmajs, bmins, bpas = [], [], []
    for stoke in beamdict.keys():
        for beamno in beamdict[stoke].keys():
            bmaj = beamdict[stoke][beamno]['beams']['BMAJarcsec']
            bmin = beamdict[stoke][beamno]['beams']['BMINarcsec']
            bpa = beamdict[stoke][beamno]['beams']['BPAdeg']

            bmajs.append(bmaj)
            bmins.append(bmin)
            bpas.append(bpa)

            beams = Beams(major=bmaj*u.arcsec, minor=bmin *
                          u.arcsec, pa=bpa*u.deg)
            beamdict[stoke][beamno].update(
                {
                    'PSFs': beams
                }
            )

    bmajs = np.array(bmajs).ravel()
    bmins = np.array(bmins).ravel()
    bpas = np.array(bpas).ravel()

    idx = (bmajs == 0) | (bmins == 0)
    bmajs[idx] = np.nan
    bmins[idx] = np.nan
    bpas[idx] = np.nan

    # Check for cutout
    cutoff = args.cutoff
    if cutoff is not None:
        cut_idx = bmajs > cutoff
        if verbose:
            print(
                f'Cutoff will mask {sum(cut_idx)} planes or {sum(cut_idx)/sum(~idx)*100}%')
        bmajs[cut_idx] = np.nan
        bmins[cut_idx] = np.nan
        bpas[cut_idx] = np.nan

    # Make mask
    mask = (np.isfinite(bmajs)) & (np.isfinite(bmins)) & (np.isfinite(bpas))

    allbeams = Beams(
        major=bmajs*u.arcsec,
        minor=bmins*u.arcsec,
        pa=bpas*u.deg
    )

    # Parse args
    target_bmaj = args.bmaj
    target_bmin = args.bmin
    target_bpa = args.bpa

    targs = [target_bmaj, target_bmin, target_bpa]

    if any([targ is None for targ in targs]) and not all([targ is None for targ in targs]):
        raise Exception('Please specify _all_ target beam params!')

    # Get common beam
    if verbose:
        print('Computing common beam...this might take some time...')
    try:
        common_beam = allbeams[mask].common_beam(tolerance=args.tolerance,
                                                 nsamps=args.nsamps,
                                                 epsilon=args.epsilon)
    except BeamError:
        if verbose:
            print("Couldn't find common beam with defaults")
            print("Trying again with smaller tolerance")

        common_beam = allbeams[mask].common_beam(tolerance=args.tolerance*0.1,
                                                 nsamps=args.nsamps,
                                                 epsilon=args.epsilon)

    if verbose:
        print('Common beam is', common_beam)

    # Set target beam
    if all([targ is None for targ in targs]):
        targbool = False
        if verbose:
            print('No target beam specified. Using smallest common beam.')
        target_beam = common_beam

    else:
        targbool = True
        target_beam = Beam(
            major=target_bmaj*u.arcsec,
            minor=target_bmin*u.arcsec,
            pa=target_bpa*u.deg
        )

    if verbose:
        print('Target beam is', target_beam)

    # Check that taget beam will deconvolve. Mask all beams that don't
    if targbool:
        if verbose:
            print('Checking that target beam will deconvolve all beams')
            print("Beams that can't will be masked")
        mask_count = 0
        for i, (beam, msk) in enumerate(
            tqdm(
                zip(allbeams, mask),
                total=len(allbeams),
                desc='Deconvolving',
                disable=(not verbose)
            )
        ):
            if not msk:
                continue
            else:
                try:
                    target_beam.deconvolve(beam)
                except ValueError:
                    mask_count += 1
                    mask[i] = False

        if verbose:
            print(
                f"{mask_count} channels or {mask_count/len(allbeams)*100}% were maksed")

    # Now comput convolving beams and masks for each beam
    if verbose:
        print('Computing and saving convolving beams for each beam')
    for stoke in tqdm(beamdict.keys(),
                      desc='Stokes',
                      disable=(not verbose)):
        for beamno in tqdm(beamdict[stoke].keys(),
                           desc='Beam',
                           disable=(not verbose)):
            beams = beamdict[stoke][beamno]['PSFs']
            if cutoff is not None:
                mask = (beams > 0) | (beams.major > cutoff*u.arcsec)
            else:
                mask = (beams > 0)
            conbeams = []
            for i, (beam, msk) in enumerate(
                zip(beams, mask)
            ):
                try:
                    conbeam = target_beam.deconvolve(beam)
                    conbeams.append(conbeam)
                except ValueError:
                    mask[i] = False
                    conbeams.append(np.nan)
            beamdict[stoke][beamno].update(
                {
                    'conbeams': conbeams,
                    'mask': mask
                }
            )

    cutdir = args.cutdir
    if cutdir is not None:
        if cutdir[-1] == '/':
            cutdir = cutdir[:-1]

    files = sorted(glob(f"{cutdir}/*/[!sm.]*.image.*.fits"))
    if files == []:
        raise Exception('No files found!')

    worker_partial = partial(
        worker,
        beamdict=beamdict,
        target_beam=target_beam,
        dryrun=args.dryrun,
        verbose=args.verbose_worker)

    list(
        tqdm(
            pool.imap(worker_partial, files),
            desc='Smoothing cutouts',
            total=len(files),
            disable=(not verbose)
        )
    )

    if verbose:
        print('Done!')


def cli():
    """Command-line interface
    """
    import argparse

    # Help string to be shown using the -h option
    descStr = """
    Smooth a field of 3D cubes to a common resolution.

    Names of output files are 'infile'.sm.fits

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        'cutdir',
        metavar='cutdir',
        type=str,
        help='Cutout directory',
    )

    parser.add_argument(
        'logdir',
        metavar='logdir',
        type=str,
        default=None,
        help='Directory containing beamlog files')

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="verbose output [False]."
    )

    parser.add_argument(
        "-vw",
        "--verbose_worker",
        dest="verbose_worker",
        action="store_true",
        help="verbose output [False].")

    parser.add_argument("-d", "--dryrun", dest="dryrun", action="store_true",
                        help="Compute common beam and stop [False].")

    parser.add_argument(
        "--bmaj",
        dest="bmaj",
        type=float,
        default=None,
        help="BMAJ to convolve to [max BMAJ from given image(s)].")

    parser.add_argument(
        "--bmin",
        dest="bmin",
        type=float,
        default=None,
        help="BMIN to convolve to [max BMAJ from given image(s)].")

    parser.add_argument(
        "--bpa",
        dest="bpa",
        type=float,
        default=None,
        help="BPA to convolve to [0].")

    parser.add_argument(
        '-m',
        '--mask',
        dest='masklist',
        type=str,
        default=None,
        help='List of channels to be masked [None]')

    parser.add_argument(
        '-c',
        '--cutoff',
        dest='cutoff',
        type=float,
        default=None,
        help='Cutoff BMAJ value -- Blank channels with BMAJ larger than this [None -- no limit]')

    parser.add_argument(
        "-t",
        "--tolerance",
        dest="tolerance",
        type=float,
        default=0.0001,
        help="tolerance for radio_beam.commonbeam.")

    parser.add_argument(
        "-e",
        "--epsilon",
        dest="epsilon",
        type=float,
        default=0.0005,
        help="epsilon for radio_beam.commonbeam.")

    parser.add_argument(
        "-n",
        "--nsamps",
        dest="nsamps",
        type=int,
        default=200,
        help="nsamps for radio_beam.commonbeam.")

    group = parser.add_mutually_exclusive_group()

    group.add_argument("--ncores", dest="n_cores", default=1,
                       type=int, help="Number of processes (uses multiprocessing).")
    group.add_argument("--mpi", dest="mpi", default=False,
                       action="store_true", help="Run with MPI.")

    args = parser.parse_args()

    verbose = args.verbose

    pool = schwimmbad.choose_pool(
        mpi=args.mpi, processes=args.n_cores)
    if args.mpi:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)


    # make it so we can use imap in serial and mpi mode
    if not isinstance(pool, schwimmbad.MultiPool):
        pool.imap = pool.map

    main(pool, args, verbose=verbose)
    pool.close()


if __name__ == "__main__":
    cli()
