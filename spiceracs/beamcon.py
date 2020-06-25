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
import psutil
print = functools.partial(print, f'[{psutil.Process().cpu_num()}]', flush=True)
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


# def getbeams(beamlog, verbose=False):
#    """
#
#    colnames=['Channel', 'BMAJarcsec', 'BMINarcsec', 'BPAdeg']
#    """
#    if verbose:
#        print(f'Getting beams from {beamlog}')
#    beams = np.genfromtxt(beamlog, names=True)
#    nchan = beams.shape[0]
#
#    return beams, nchan


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


def getfacs(datadict, new_beam, verbose=False):
    """Get beam info
    """
    conbms = []
    facs = []
    for oldbeam in datadict['oldbeams']:
        if verbose:
            print(f"Current beam is", oldbeam)
        conbm = new_beam.deconvolve(oldbeam)
        fac, amp, outbmaj, outbmin, outbpa = au2.gauss_factor(
            [
                conbm.major.to(u.arcsec).value,
                conbm.minor.to(u.arcsec).value,
                conbm.pa.to(u.deg).value
            ],
            beamOrig=[
                oldbeam.major.to(u.arcsec).value,
                oldbeam.minor.to(u.arcsec).value,
                oldbeam.pa.to(u.deg).value
            ],
            dx1=datadict['dx'].to(u.arcsec).value,
            dy1=datadict['dy'].to(u.arcsec).value
        )
        conbms.append(conbm)
        facs.append(fac)
    return conbms, facs


def smooth(image, dy, conbeam, sfactor, verbose=False):
    """Do the smoothing
    """
    if np.isnan(conbeam):
        return image*np.nan
    if np.isnan(image).all():
        return image
    else:
        # using Beams package
        if verbose:
            print(f'Using convolving beam', conbeam)
        pix_scale = dy
        gauss_kern = conbeam.as_kernel(dy)

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


def worker(idx, start, cubedict):
    cube = SpectralCube.read(cubedict["filename"])
    plane = cube.unmasked_data[start+idx].value
    newim = smooth(plane, cubedict['dy'], cubedict['conbeams']
                   [start+idx], cubedict['sfactors'][start+idx], verbose=False)
    return newim


def main(pool, args, verbose=True):
    if args.dryrun:
        if verbose:
            print('Doing a dry run -- no files will be saved')
    # Fix up outdir
    if args.mpi:
        from mpi4py import MPI
        mpiComm = MPI.COMM_WORLD
        n_cores = mpiComm.Get_size()
        #mpiRank = mpiComm.Get_rank()

    logdir = args.logdir
    if logdir is not None:
        if logdir[-1] == '/':
            logdir = logdir[:-1]

    beamdict = getbeams(logdir)

    bmajs, bmins, bpas = [], [], []
    for stoke in beamdict.keys():
        for beamno in beamdict[stoke].keys():
            bmaj = beamdict[stoke][beamno]['beams']['BMAJarcsec']
            bmin = beamdict[stoke][beamno]['beams']['BMINarcsec']
            bpa = beamdict[stoke][beamno]['beams']['BPAdeg']

            bmajs.append(bmaj)
            bmins.append(bmin)
            bpas.append(bpa)

            beams = Beams(major=bmaj*u.arcsec, minor=bmin*u.arcsec, pa=bpa*u.deg)
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
    bpas[idx]  = np.nan

    cutoff = args.cutoff
    if cutoff is not None:
        cut_idx = bmajs > cutoff
        if verbose:
            print(f'Cutoff will mask {sum(cut_idx)} planes or {sum(cut_idx)/sum(~idx)*100}%')
        bmajs[cut_idx] = np.nan
        bmins[cut_idx] = np.nan
        bpas[cut_idx]  = np.nan
    
    mask = (np.isfinite(bmajs)) & (np.isfinite(bmins)) & (np.isfinite(bpas))

    allbeams = Beams(
        major = bmajs*u.arcsec,
        minor = bmins*u.arcsec,
        pa = bpas*u.deg
    )

    # Parse args
    target_bmaj = args.bmaj
    target_bmin = args.bmin
    target_bpa = args.bpa

    targs = [target_bmaj, target_bmin, target_bpa]

    if any([targ is None for targ in targs]) and not all([targ is None for targ in targs]):
        raise Exception('Please specify _all_ target beam params!')

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

    if all([targ is None for targ in targs]):
        targbool = False
        if verbose:
            print('No target beam specified. Using smallest common beam.')
        target_beam = common_beam

    else:
        targbool = True
        target_beam = Beam(
            major = target_bmaj*u.arcsec,
            minor = target_bmin*u.arcsec,
            pa = target_bpa*u.deg
        )

    if verbose:
        print('Target beam is', target_beam)

    embed()
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
            print(f"{mask_count} channels or {mask_count/len(allbeams)*100}% were maksed")


    cutdir = args.cutdir
    if cutdir is not None:
        if cutdir[-1] == '/':
            cutdir = cutdir[:-1]

    files = sorted(glob(f"{cutdir}/*/*.image.*.fits"))
    if files == []:
        raise Exception('No files found!')

    
    beams = []
    nchans = []
    datadict = {}
    masks = []
    for i, file in enumerate(tqdm(files)):
        # Set up files
        datadict[f"cube_{i}"] = {}
        datadict[f"cube_{i}"]["filename"] = file
        # Get metadata
        with fits.open(file, memmap=True, mode='denywrite') as hdulist:
            header = hdulist[0].header
            dxas = header['CDELT1']*-1*u.deg
            datadict[f"cube_{i}"]["dx"] = dxas
            dyas = header['CDELT2']*u.deg
            datadict[f"cube_{i}"]["dy"] = dyas
            # Get beam info
            dirname = os.path.dirname(file)
            datadict[f"cube_{i}"]["outdir"] = dirname
            filename = os.path.basename(file)
            basename = filename[filename.find('image.'):]
            if dirname == '':
                dirname = '.'
            #beamlog = f"{logdir}/beamlog.{basename}".replace('.fits', '.txt')
            stoke = basename[15+basename.find('image.restored.')
                                              :15+basename.find('image.restored.')+1]
            beamno = basename[basename.find(
                'beam')+4:basename.find('beam')+4+2]
            datadict[f"cube_{i}"]["beamlog"] = beamdict[stoke][beamno]['beamlog']

            beam = beamdict[stoke][beamno]['beams']
            nchan = beamdict[stoke][beamno]['nchan']

            if args.cutoff is not None:
                frac = len([beam['BMAJarcsec'].ravel() > args.cutoff]) / \
                    len(beam['BMAJarcsec'].ravel())
                # if verbose:
                # print(
                #    f'Cutoff will blank {round(frac*100,2)}% of channels in {basename}')
            # Find bad chans
            #cube = SpectralCube.read(file)
            #mask = cube[:, cube.shape[1]//2, cube.shape[2]//2].mask.view()
            fitscube = hdulist[0].data
            subcube = fitscube[:, 0, fitscube.shape[2] //
                               2, fitscube.shape[3]//2]
            mask = (subcube == 0) | (np.isnan(subcube))
            masks.append(mask)
            # Record beams
            beams.append(beam)
            nchans.append(nchan)
    beams = np.array(beams)
    nchans = np.array(nchans)
    # Do dome masking
    beams['BMAJarcsec'][beams['BMAJarcsec'] == 0] = np.nan
    beams['BMINarcsec'][beams['BMAJarcsec'] == 0] = np.nan
    beams['BMINarcsec'][beams['BMINarcsec'] == 0] = np.nan
    beams['BMAJarcsec'][beams['BMINarcsec'] == 0] = np.nan

    if args.cutoff is not None:
        frac = len([beams['BMAJarcsec'].ravel() > args.cutoff]) / \
            len(beams['BMAJarcsec'].ravel())
        if verbose:
            print(
                f'Cutoff will blank {round(frac*100,2)}% of channels, in all cubes')
        beams['BMAJarcsec'][beams['BMAJarcsec'] > args.cutoff] = np.nan
        beams['BMINarcsec'][beams['BMAJarcsec'] > args.cutoff] = np.nan

    totalmask = sum(masks) > 0

    for i, _ in enumerate(beams['BMAJarcsec']):
        beams['BMAJarcsec'][i][totalmask] = np.nan
        beams['BMINarcsec'][i][totalmask] = np.nan
        beams['BPAdeg'][i][totalmask] = np.nan
        datadict[f"cube_{i}"]["oldbeams"] = Beams(beams['BMAJarcsec'][i].ravel(
        )*u.arcsec, beams['BMINarcsec'][i].ravel()*u.arcsec, beams['BPAdeg'][i].ravel()*u.deg)

    if args.masklist is not None:
        masklist = np.loadtxt(args.masklist) == 1
        for i, _ in enumerate(beams['BMAJarcsec']):
            beams['BMAJarcsec'][i][masklist] = np.nan
            beams['BMINarcsec'][i][masklist] = np.nan
            beams['BPAdeg'][i][masklist] = np.nan
            datadict[f"cube_{i}"]["oldbeams"] = Beams(beams['BMAJarcsec'][i].ravel(
            )*u.arcsec, beams['BMINarcsec'][i].ravel()*u.arcsec, beams['BPAdeg'][i].ravel()*u.deg)

        # for chan in masklist:

    if not all(elem == nchans[0] for elem in nchans):
        raise Exception('Unequal channel count in beamlogs!')

    beamlst = Beams(beams['BMAJarcsec'].ravel(
    )*u.arcsec, beams['BMINarcsec'].ravel()*u.arcsec, beams['BPAdeg'].ravel()*u.deg)

    try:
        big_beam = beamlst[~np.isnan(beamlst)].common_beam(tolerance=args.tolerance,
                                                           nsamps=args.nsamps,
                                                           epsilon=args.epsilon)
    except BeamError:
        if verbose:
            print("Couldn't find common beam with defaults")
            print("Trying again with smaller tolerance")

        big_beam = beamlst[~np.isnan(beamlst)].common_beam(tolerance=args.tolerance*0.1,
                                                           nsamps=args.nsamps,
                                                           epsilon=args.epsilon)

    if verbose:
        print(f'Smallest common beam is', big_beam)
    # Parse args
    bmaj = args.bmaj
    bmin = args.bmin
    bpa = args.bpa

    # Set to largest
    if bpa is None and bmin is None and bmaj is None:
        bpa = big_beam.pa.to(u.deg)
    else:
        bpa = 0*u.deg
    if bmaj is None:
        bmaj = round_up(big_beam.major.to(u.arcsec))
    elif bmaj*u.arcsec < round_up(big_beam.major.to(u.arcsec)):
        raise Exception('Selected BMAJ is too small!')
    else:
        bmaj *= u.arcsec
    if bmin is None:
        bmin = round_up(big_beam.minor.to(u.arcsec))
    elif bmin*u.arcsec < round_up(big_beam.minor.to(u.arcsec)):
        raise Exception('Selected BMIN is too small!')
    else:
        bmin *= u.arcsec

    new_beam = Beam(
        bmaj,
        bmin,
        bpa
    )
    if verbose:
        print(f'Final beam is', new_beam)

    if not args.dryrun:
        for key in tqdm(datadict.keys(), desc='Working on cubes separately'):
            conbms, facs = getfacs(datadict[key], new_beam, verbose=False)
            cube = SpectralCube.read(datadict[key]["filename"])
            # Set up output file
            outname = "sm." + os.path.basename(datadict[key]["filename"])
            outdir = datadict[key]["outdir"]
            outfile = f'{outdir}/{outname}'
            if verbose:
                print(f'Initialsing to {outfile}')
            if not os.path.isfile(outfile):
                copyfile(datadict[key]["filename"], outfile, verbose=True)

            cubedict = datadict[key]
            cubedict["conbeams"] = conbms
            cubedict["sfactors"] = facs

            if not args.mpi:
                n_cores = args.n_cores
            width_max = n_cores
            width = cpu_to_use(width_max, cube.shape[0])
            n_chunks = cube.shape[0]//width

            for i in trange(
                    n_chunks, disable=(not verbose),
                    desc='Smoothing in chunks'
            ):
                start = i*width
                stop = start+width

                func = functools.partial(
                    worker, start=start, cubedict=cubedict)
                arr_out = list(pool.map(func, [idx for idx in range(width)]))
                arr_out = np.array(arr_out)

                with fits.open(outfile, mode='update', memmap=True) as outfh:
                    outfh[0].data[start:stop, 0, :, :] = arr_out[:]
                    outfh.flush()

            if verbose:
                print('Updating header...')
            with fits.open(outfile, mode='update', memmap=True) as outfh:
                outfh[0].header = new_beam.attach_to_header(outfh[0].header)
                outfh.flush()
            # print(arr_out)
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

    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
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
