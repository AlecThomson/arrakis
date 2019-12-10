#!/usr/bin/env python
from spiceracs.utils import getdata, copyfile, getfreq, cpu_to_use
from RMtools_3D.do_RMsynth_3D import run_rmsynth
from spectral_cube import SpectralCube
import warnings
import numpy as np
from astropy.io import fits
from tqdm import tqdm, trange
import os
import pdb
import pymongo
import multiprocessing as mp
import ctypes as c
import time
import functools
print = functools.partial(print, flush=True)


class moments:
    """Methods for produce moment maps.
    """

    def __init__(self, cube, basename, n_cores, verbose=True):
        self.cube = cube
        self.verbose = verbose
        self.basename = basename
        self.n_cores = n_cores

    def mu(self, outdir='.', stoke=''):
        """Mean moment - freq axis first
        """
        def mu_worker(queue, inarr, outarr):
            while not queue.empty():
                start, idx = queue.get()
                yav = np.nanmean(np.array(self.cube[:, i, :]), axis=0)
                yav[yav == 0] = np.nan
                outarr[idx, :] = yav[:]

        try:
            os.mkdir(f'{outdir}/moments/')
            print('Made directory.')
        except FileExistsError:
            print('Directory exists.')
        momfilename = self.basename.replace(
            '.i.', f'.{stoke}.').replace('contcube', 'mu')
        momfile = f'{outdir}/moments/{momfilename}'
        blank = self.cube[0]
        blank.write(momfile, overwrite=True, format='fits')

        if self.verbose:
            print(f'Writing to {momfile}...')
        procs = []
        width_max = 100
        width = cpu_to_use(width_max, self.cube.shape[1])
        n_chunks = self.cube.shape[1]//width

        outshape = list(self.cube.shape)
        outshape[1] = width

        # shared, can be used from multiple processes
        mp_arr_out = mp.Array(c.c_double, int(np.prod(outshape)))
        # then in each new process create a new numpy array using:
        # mp_arr and arr share the same memory
        buff_arr_out = np.frombuffer(mp_arr_out.get_obj())
        # make it two-dimensional
        # b and arr share the same memory
        arr_out = buff_arr_out.reshape(outshape)
        arr_out[:] = np.zeros(outshape)*np.nan

        inshape = list(self.cube.shape)
        inshape[1] = width

        mp_arr_in = mp.Array(c.c_double, int(np.prod(inshape)))
        # then in each new process create a new numpy array using:
        # mp_arr and arr share the same memory
        buff_arr_in = np.frombuffer(mp_arr_in.get_obj())
        # make it two-dimensional
        # b and arr share the same memory
        arr_in = buff_arr_in.reshape(inshape)
        arr_in[:] = np.zeros(inshape)*np.nan

        q = mp.Queue()

        tic = time.perf_counter()
        for i in trange(
                n_chunks, disable=(not self.verbose),
                desc='Making mu in chunks'
        ):
            # if verbose:
                #print('Beginning multiprocessing...')
            procs = []
            for j in range(self.n_cores):
                proc = mp.Process(target=mu_worker, args=(
                    q, arr_in, arr_out
                )
                )
                procs.append(proc)

            start = i*width
            stop = start+width
            # if verbose:
            #print('Copying chunk of data to memory...')

            dat = np.array(self.cube[:, start:stop, :])
            dat[dat == 0] = np.nan
            arr_in[:] = dat
            for idx in range(stop-start):
                q.put([start, idx])
            # if verbose:
                # print('Processing...')
            for p in procs:
                p.start()
            for p in procs:
                p.join()
            for p in procs:
                p.terminate()
            # if verbose:
                #print('Writing output to file...')
            with fits.open(momfile, mode='update', memmap=True) as outfh:
                outfh[0].data[start:stop] = arr_out[:]
                outfh.flush()

        toc = time.perf_counter()

        if self.verbose:
            print(f'Time taken was {toc - tic}s')

    def sigma(self, outdir='.', stoke=''):
        """
        Standard deviation moment - freq axis first
        """
        def sigma_worker(queue, inarr, outarr):
            while not queue.empty():
                start, idx = queue.get()
                ystd = np.nanstd(np.array(inarr[:, idx, :]), axis=0, ddof=1)
                ystd[ystd == 0] = np.nan
                outarr[idx, :] = ystd[:]

        try:
            os.mkdir(f'{outdir}/moments/')
            print('Made directory.')
        except FileExistsError:
            print('Directory exists.')

        momfilename = self.basename.replace(
            '.i.', f'.{stoke}.').replace('contcube', 'sigma')
        momfile = f'{outdir}/moments/{momfilename}'
        blank = self.cube[0]
        blank.write(momfile, overwrite=True, format='fits')

        if self.verbose:
            print(f'Writing to {momfile}...')
        procs = []
        width_max = 100
        width = cpu_to_use(width_max, self.cube.shape[1])
        n_chunks = self.cube.shape[1]//width

        outshape = list(self.cube.shape)
        outshape[1] = width

        # shared, can be used from multiple processes
        mp_arr_out = mp.Array(c.c_double, int(np.prod(outshape)))
        # then in each new process create a new numpy array using:
        # mp_arr and arr share the same memory
        buff_arr_out = np.frombuffer(mp_arr_out.get_obj())
        # make it two-dimensional
        # b and arr share the same memory
        arr_out = buff_arr_out.reshape(outshape)
        arr_out[:] = np.zeros(outshape)*np.nan

        inshape = list(self.cube.shape)
        inshape[1] = width

        mp_arr_in = mp.Array(c.c_double, int(np.prod(inshape)))
        # then in each new process create a new numpy array using:
        # mp_arr and arr share the same memory
        buff_arr_in = np.frombuffer(mp_arr_in.get_obj())
        # make it two-dimensional
        # b and arr share the same memory
        arr_in = buff_arr_in.reshape(inshape)
        arr_in[:] = np.zeros(inshape)*np.nan

        q = mp.Queue()

        tic = time.perf_counter()
        for i in trange(
                n_chunks, disable=(not self.verbose),
                desc='Making sigma in chunks'
        ):
            # if verbose:
                #print('Beginning multiprocessing...')
            procs = []
            for j in range(self.n_cores):
                proc = mp.Process(target=sigma_worker, args=(
                    q, arr_in, arr_out
                )
                )
                procs.append(proc)

            start = i*width
            stop = start+width
            # if verbose:
            #print('Copying chunk of data to memory...')

            dat = np.array(self.cube[:, start:stop, :])
            dat[dat == 0] = np.nan
            arr_in[:] = dat
            for idx in range(stop-start):
                q.put([start, idx])
            # if verbose:
                # print('Processing...')
            for p in procs:
                p.start()
            for p in procs:
                p.join()
            for p in procs:
                p.terminate()
            # if verbose:
                #print('Writing output to file...')
            with fits.open(momfile, mode='update', memmap=True) as outfh:
                outfh[0].data[start:stop] = arr_out[:]
                outfh.flush()

        toc = time.perf_counter()

        if self.verbose:
            print(f'Time taken was {toc - tic}s')


def momentloop(datadict, n_cores, outdir='.', verbose=True):
    """Produce moment maps.

    Args:
        datadict (dict): Dictionary containing tables and cubes.

    Kwargs:
        outdir (str): Directory to save data to.
        verbose (bool): Print out messages.

    Returns:

    """
    basename = os.path.basename(datadict['i_file'])
    for stokes in ['p', 'q', 'u']:
        if verbose:
            print(f'Making moments for Stokes {stokes}...')
        mom = moments(datadict[f'{stokes}_cube'],
                      basename, n_cores, verbose=verbose)
        if verbose:
            print(f'Making mean...')
        mom.mu(outdir=outdir, stoke=stokes)
        if verbose:
            print(f'Making std...')
        mom.sigma(outdir=outdir, stoke=stokes)


def makepi_worker(queue, inarr_q, inarr_u, outarr, verbose):
    while not queue.empty():
        start, idx = queue.get()
        ypi = np.hypot(
            np.expand_dims(inarr_q[:, idx, :], axis=1),
            np.expand_dims(inarr_u[:, idx, :], axis=1)
        )
        ypi[ypi == 0] = np.nan
        outarr[:, :, idx, :] = ypi[:]


def makepi(datadict, n_cores, outdir='.', verbose=True):
    """Make a polarized itensity cube.

    Args:
        datadict (dict): Dictionary containing tables and cubes.

    Kwargs:
        verbose (bool): Print out messages.

    Returns:
        datadict (dict): Dictionary containing tables and cubes, now
        updated with PI data.

    """
    if outdir[-1] == '/':
        outdir = outdir[:-1]
    try:
        os.mkdir(f'{outdir}/moments/')
        print('Made directory.')
    except FileExistsError:
        print('Directory exists.')

    pifilename = os.path.basename(datadict['i_file'].replace('.i.', '.p.'))
    pifile = f'{outdir}/moments/{pifilename}'

    try:
        p_cube = SpectralCube.read(pifile, mode='denywrite')
        datadict['p_file'] = pifile
        datadict['p_cube'] = p_cube
        return datadict

    except FileNotFoundError:
        if verbose:
            print(f'Writing to {pifile}...')
        try:
            with open(pifile, 'rb') as fh:
                pass
        except FileNotFoundError:
            copyfile(datadict['i_file'], pifile, verbose=verbose)

        if verbose:
            print('Initialising arrays...')
            print('Note: Output data will be held in memory.')

        procs = []
        width_max = 100
        width = cpu_to_use(width_max, datadict['i_cube'].shape[1])
        n_chunks = datadict['i_cube'].shape[1]//width

        outshape = list(datadict['i_cube'].shape)
        outshape.insert(1, 1)  # For Stokes axis
        outshape[2] = width

        # shared, can be used from multiple processes
        mp_arr_out = mp.Array(c.c_double, int(np.prod(outshape)))
        # then in each new process create a new numpy array using:
        # mp_arr and arr share the same memory
        buff_arr_out = np.frombuffer(mp_arr_out.get_obj())
        # make it two-dimensional
        # b and arr share the same memory
        arr_out = buff_arr_out.reshape(outshape)
        arr_out[:] = np.zeros(outshape)*np.nan

        inshape = list(datadict['i_cube'].shape)
        inshape[1] = width

        mp_arr_in_q = mp.Array(c.c_double, int(np.prod(inshape)))
        # then in each new process create a new numpy array using:
        # mp_arr and arr share the same memory
        buff_arr_in_q = np.frombuffer(mp_arr_in_q.get_obj())
        # make it two-dimensional
        # b and arr share the same memory
        arr_in_q = buff_arr_in_q.reshape(inshape)
        arr_in_q[:] = np.zeros(inshape)*np.nan

        mp_arr_in_u = mp.Array(c.c_double, int(np.prod(inshape)))
        # then in each new process create a new numpy array using:
        # mp_arr and arr share the same memory
        buff_arr_in_u = np.frombuffer(mp_arr_in_u.get_obj())
        # make it two-dimensional
        # b and arr share the same memory
        arr_in_u = buff_arr_in_u.reshape(inshape)
        arr_in_u[:] = np.zeros(inshape)*np.nan

        q = mp.Queue()

        tic = time.perf_counter()
        for i in trange(
                n_chunks, disable=(not verbose), desc='Making PI cube in chunks'
        ):
            # if verbose:
                #print('Beginning multiprocessing...')
            procs = []
            for j in range(n_cores):
                proc = mp.Process(target=makepi_worker, args=(
                    q, arr_in_q, arr_in_u, arr_out, verbose
                )
                )
                procs.append(proc)

            start = i*width
            stop = start+width
            # if verbose:
            #print('Copying chunk of data to memory...')

            datq = np.array(datadict['q_cube'][:, start:stop, :])
            datu = np.array(datadict['u_cube'][:, start:stop, :])
            datq[datq == 0] = np.nan
            datu[datu == 0] = np.nan
            arr_in_q[:] = datq
            arr_in_u[:] = datu
            for idx in range(stop-start):
                q.put([start, idx])
            # if verbose:
                # print('Processing...')
            for p in procs:
                p.start()
            #for p in procs:
            #    p.join()
            #for p in procs:
            #    p.terminate()
            # if verbose:
                #print('Writing output to file...')
            with fits.open(pifile, mode='update', memmap=True) as outfh:
                outfh[0].data[:, :, start:stop] = arr_out[:]
                outfh.flush()

        toc = time.perf_counter()

        if verbose:
            print('PI cube done!')
            print(f'Time taken was {toc - tic}s')

        p_cube = SpectralCube.read(pifile, mode='denywrite')
        datadict['p_file'] = pifile
        datadict['p_cube'] = p_cube

        return datadict


def makezero_worker(queue, inarr_q, inarr_u, inarr_f, outarr, verbose):
    while not queue.empty():
        start, idx = queue.get()
        dataArr = run_rmsynth(
            np.array(inarr_q[:, idx, :]),
            np.array(inarr_q[:, idx, :]),
            np.array(inarr_f), phiMax_radm2=1000, nSamples=5.0,
            weightType="uniform", fitRMSF=False, nBits=32,
            verbose=False, not_rmsf=True
        )
        FDFcube, phiArr_radm2, lam0Sq_m2, lambdaSqArr_m2 = dataArr
        dphi = np.diff(phiArr_radm2)[0]
        mom0 = np.nansum(abs(FDFcube)*dphi, axis=0)
        mom0[mom0 == 0] = np.nan
        outarr[start+idx] = mom0[:]


def makezero(datadict, n_cores, outdir='.', verbose=True):
    """Make a Faraday moment 0 cube.

    Args:
        datadict (dict): Dictionary containing tables and cubes.

    Kwargs:
        verbose (bool): Print out messages.

    Returns:
        datadict (dict): Dictionary containing tables and cubes, now
        updated with PI data.

    """
    if outdir[-1] == '/':
        outdir = outdir[:-1]

    if verbose:
        print(f'Tring to make {outdir}/moments/...')

    try:
        os.mkdir(f'{outdir}/moments/')
        print('Made directory.')
    except FileExistsError:
        print('Directory exists.')

    momfilename = os.path.basename(datadict['i_file'].replace(
        '.i.', f'.p.').replace('contcube', 'mom0'))
    momfile = f'{outdir}/moments/{momfilename}'

    blank = datadict['i_cube'][0]
    if verbose:
        print(f'Writing to {momfile}...')
    blank.write(momfile, overwrite=True, format='fits')

    if verbose:
        print('Initialising arrays...')
        print('Note: Output data will be held in memory.')

    freq = np.array(getfreq(datadict['q_cube']))
    freq_arr = mp.Array(c.c_double, len(freq))
    freq_arr[:] = freq[:]

    outshape = datadict['q_cube'].shape[1:]
    # shared, can be used from multiple processes
    mp_arr_out = mp.Array(c.c_double, int(np.prod(outshape)))
    # then in each new process create a new numpy array using:
    # mp_arr and arr share the same memory
    buff_arr_out = np.frombuffer(mp_arr_out.get_obj())
    # make it two-dimensional
    arr_out = buff_arr_out.reshape(outshape)  # b and arr share the same memory
    arr_out[:] = np.zeros(outshape)*np.nan

    procs = []
    width_max = 100
    width = cpu_to_use(width_max, datadict['q_cube'].shape[1])
    n_chunks = datadict['q_cube'].shape[1]//width

    inshape = list(datadict['q_cube'].shape)
    inshape[1] = width

    mp_arr_in_q = mp.Array(c.c_double, int(np.prod(inshape)))
    # then in each new process create a new numpy array using:
    # mp_arr and arr share the same memory
    buff_arr_in_q = np.frombuffer(mp_arr_in_q.get_obj())
    # make it two-dimensional
    # b and arr share the same memory
    arr_in_q = buff_arr_in_q.reshape(inshape)
    arr_in_q[:] = np.zeros(inshape)*np.nan

    mp_arr_in_u = mp.Array(c.c_double, int(np.prod(inshape)))
    # then in each new process create a new numpy array using:
    # mp_arr and arr share the same memory
    buff_arr_in_u = np.frombuffer(mp_arr_in_u.get_obj())
    # make it two-dimensional
    # b and arr share the same memory
    arr_in_u = buff_arr_in_u.reshape(inshape)
    arr_in_u[:] = np.zeros(inshape)*np.nan

    q = mp.Queue()

    tic = time.perf_counter()
    for i in trange(
            n_chunks, disable=(not verbose), desc='Making M0 in chunks'
    ):
        # if verbose:
            #print('Beginning multiprocessing...')
        procs = []
        for j in range(n_cores):
            proc = mp.Process(target=makezero_worker, args=(
                q, arr_in_q, arr_in_u, freq_arr, arr_out, verbose
            )
            )
            procs.append(proc)

        start = i*width
        stop = start+width
        # if verbose:
        #print('Copying chunk of data to memory...')
        datq = np.array(datadict['q_cube'][:, start:stop, :])
        datu = np.array(datadict['u_cube'][:, start:stop, :])
        datq[datq == 0] = np.nan
        datu[datu == 0] = np.nan
        arr_in_q[:] = datq
        arr_in_u[:] = datu
        for idx in range(stop-start):
            q.put([start, idx])
        # if verbose:
            # print('Processing...')
        for p in procs:
            p.start()
        for p in procs:
            p.join()
        for p in procs:
            p.terminate()
        # if verbose:
            #print('Writing output to file...')
        with fits.open(momfile, mode='update', memmap=True) as outfh:
            outfh[0].data[start:stop] = arr_out[start:stop]
            outfh.flush()
    toc = time.perf_counter()
    if verbose:
        print('Moment 0 done!')
        print(f'Time taken was {toc - tic}s')


def main(args, verbose=True):
    """Main script.
    """
    # Sort out args
    cubedir = args.cubedir
    tabledir = args.tabledir
    mapdir = args.mapdir
    outdir = args.outdir

    n_cores = args.n_cores
    if args.n_cores is None:
        n_cores = mp.cpu_count()

    # Read in data
    if verbose:
        print('Reading data...')
    datadict = getdata(cubedir, tabledir, mapdir, verbose=verbose)

    # Make PI cube
    if args.picube:
        if verbose:
            print('Making polarized intensity cube...')
        datadict = makepi(datadict, n_cores, outdir=outdir, verbose=verbose)

    # Compute moments
    if args.farnes:
        if verbose:
            print('Making moment maps...')
        momentloop(datadict, n_cores, outdir=outdir, verbose=verbose)

    if args.zero:
        if verbose:
            print('Making Faraday zeroth moment...')
        makezero(datadict, n_cores, outdir=outdir, verbose=verbose)

    if verbose:
        print('Done!')


def cli():
    """Command-line interface.
    """
    import argparse
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
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
    SPICE-RACS Stage 3:
    Produce 'Faraday' moments. Files will be saved to 'outdir/moments/'.

    Note: A PI cube will also be produced.

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
        'mapdir',
        metavar='mapdir',
        type=str,
        help='Directory containing image data in FITS format.')

    parser.add_argument(
        'outdir',
        metavar='outdir',
        type=str,
        help='Directory to store moment maps.')

    parser.add_argument(
        "-v",
        dest="verbose",
        action="store_true",
        help="Verbose output [False]."
    )

    parser.add_argument(
        "--picube",
        dest="picube",
        action="store_true",
        help="Make a PI cube [False]."
    )

    parser.add_argument(
        "--farnes",
        dest="farnes",
        action="store_true",
        help="Make Farnes (2018) moments [False]."
    )

    parser.add_argument(
        "--zero",
        dest="zero",
        action="store_true",
        help="Make Zeroth moment using RM synthesis [False]."
    )

    parser.add_argument(
        "--ncores",
        dest="n_cores",
        type=int, help="Number of processes [use all available].")

    args = parser.parse_args()
    verbose = args.verbose

    main(args, verbose=verbose)


if __name__ == "__main__":
    cli()
