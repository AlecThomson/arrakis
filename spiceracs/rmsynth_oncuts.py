#!/usr/bin/env python
from spiceracs.utils import getfreq
import numpy as np
import os
import pymongo
import sys
import subprocess
import time
from tqdm import tqdm, trange
import warnings
from RMtools_1D import do_RMsynth_1D
from RMtools_3D import do_RMsynth_3D
from spectral_cube import SpectralCube
from astropy.io import fits
import matplotlib.pyplot as plt
from RMutils.util_misc import create_frac_spectra


def rmsythoncut(args):
    i, clargs, outdir, freqfile, verbose = args

    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

    # Basic querey
    if clargs.pol and not clargs.unres:
        myquery = {"polarized": True}
    if clargs.unres and not clargs.pol:
        myquery = {"resolved": False}
    if clargs.pol and clargs.unres:
        myquery = {"$and": [{"resolved": False}, {"polarized": True}]}
    else:
        myquery = {}

    doc = mycol.find(myquery).sort("flux_peak", -1)

    iname = doc[i]['island_name']
    qfile = f"{outdir}/{doc[i]['q_file']}"
    ufile = f"{outdir}/{doc[i]['u_file']}"
    vfile = f"{outdir}/{doc[i]['v_file']}"

    with fits.open(vfile) as hdulist:
        rms = hdulist[0].data


    if clargs.prefixOut is not None:
        prefix = clargs.prefixOut + iname
    else:
        prefix = iname

    header, dataQ = do_RMsynth_3D.readFitsCube(qfile, clargs.rm_verobse)

    rmsArr = np.ones_like(dataQ)*rms

    # Run RM-synthesis on the cubes
    dataArr = do_RMsynth_3D.run_rmsynth(
        dataQ=dataQ,
        dataU=do_RMsynth_3D.readFitsCube(ufile, clargs.rm_verobse)[1],
        freqArr_Hz=do_RMsynth_3D.readFreqFile(freqfile, clargs.rm_verobse),
        dataI=None,
        rmsArr=rmsArr,
        phiMax_radm2=clargs.phiMax_radm2,
        dPhi_radm2=clargs.dPhi_radm2,
        nSamples=clargs.nSamples,
        weightType=clargs.weightType,
        fitRMSF=clargs.fitRMSF,
        nBits=32,
        verbose=clargs.rm_verobse,
        not_rmsf=clargs.not_RMSF)

    # Write to files
    do_RMsynth_3D.writefits(dataArr,
                            headtemplate=header,
                            fitRMSF=False,
                            prefixOut=prefix,
                            outDir=outdir,
                            write_seperate_FDF=clargs.write_seperate_FDF,
                            not_rmsf=clargs.not_RMSF,
                            nBits=32,
                            verbose=verbose)

    myquery = {"island_name": iname}
    newvalues = {"$set": {"rmsynth": True}}
    mycol.update_one(myquery, newvalues)


def rmsythoncut_i(args):
    i, clargs, freq, outdir, verbose = args

    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

    # Basic querey
    if clargs.pol and not clargs.unres:
        myquery = {"polarized": True}
    if clargs.unres and not clargs.pol:
        myquery = {"resolved": False}
    if clargs.pol and clargs.unres:
        myquery = {"$and": [{"resolved": False}, {"polarized": True}]}
    else:
        myquery = {}

    doc = mycol.find(myquery).sort("flux_peak", -1)

    iname = f"{outdir}/{doc[i]['island_name']}.validate"
    data = np.array(SpectralCube.read(doc[i]['i_file']))
    mom = np.nansum(data, axis=0)
    idx = np.unravel_index(np.argmax(mom), mom.shape)

    plt.ion()
    plt.figure()
    plt.imshow(mom, origin='lower', cmap='cubehelix_r')
    plt.scatter(idx[1], idx[0], c='r', marker='x')
    plt.show()
    _ = input("Press [enter] to continue")  # wait for input from the user
    plt.close()    # close the figure to show the next one.

    data = data[:, idx[0], idx[1]]

    with fits.open(doc[i]['v_file']) as hdulist:
        noise = hdulist[0].data

    plt.ion()
    plt.figure()
    plt.step(freq/1e9, data)

    imod, qArr, uArr, dqArr, duArr, fitDict = \
        create_frac_spectra(freqArr=freq,
                            IArr=data,
                            QArr=data,
                            UArr=data,
                            dIArr=noise,
                            dQArr=noise,
                            dUArr=noise,
                            polyOrd=3,
                            verbose=True,
                            debug=False)
    plt.plot(freq/1e9, imod)
    plt.xlabel(r'$\nu$ [GHz]')
    plt.ylabel(r'$I$ [Jy/beam]')
    plt.tight_layout()
    plt.show()
    _ = input("Press [enter] to continue")  # wait for input from the user
    plt.close()    # close the figure to show the next one.

    data = data - imod
    data = data - np.nanmean(data)
    plt.ion()
    plt.figure()
    plt.step(freq/1e9, data)
    plt.xlabel(r'$\nu$ [GHz]')
    plt.ylabel(r'$I-\mathrm{model}(I)-\mathrm{mean}(\mathrm{model}(I))$')
    plt.tight_layout()
    plt.show()
    _ = input("Press [enter] to continue")  # wait for input from the user
    plt.close()    # close the figure to show the next one.

    datalist = [freq, data, data, dqArr, duArr]

    nSamples = clargs.nSamples
    phi_max = clargs.phiMax_radm2
    mDict, aDict = do_RMsynth_1D.run_rmsynth(
        datalist,
        phiMax_radm2=phi_max,
        nSamples=nSamples,
        verbose=True
    )
    plt.ion()
    plt.figure()
    plt.plot(aDict['phiArr_radm2'], abs(aDict['dirtyFDF']))
    plt.xlabel(r'$\phi$ [rad m$^{-2}$]')
    plt.ylabel(r'Dirty FDF (Stokes I)')
    plt.tight_layout()
    plt.show()
    _ = input("Press [enter] to continue")  # wait for input from the user
    plt.close()    # close the figure to show the next one.
    do_RMsynth_1D.saveOutput(mDict, aDict, iname, verbose=verbose)


def main(pool, args, verbose=False):
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
    if args.unres and not args.pol:
        myquery = {"resolved": False}
    if args.pol and args.unres:
        myquery = {"$and": [{"resolved": False}, {"polarized": True}]}
    else:
        myquery = {}

    mydoc = mycol.find(myquery).sort("flux_peak", -1)
    count = mycol.count_documents(myquery)

    if args.limit is not None:
        count = args.limit

    # Make frequency file
    freq, freqfile = getfreq(mydoc[0]['q_file'], outdir=os.path.dirname(
        mydoc[0]['q_file']), filename='frequencies.txt', verbose=verbose)
    freq = np.array(freq)
    if verbose:
        print(f'Running RMsynth on {count} sources')

    if args.validate:
        inputs = [[i, args, freq, outdir, verbose] for i in range(count)]
        if (pool.__class__.__name__ is 'MPIPool' or
                pool.__class__.__name__ is 'SerialPool'):
            if verbose:
                print('Running RM synth on Stokes I...')
            tic = time.perf_counter()
            list(pool.map(rmsythoncut_i, inputs))
            toc = time.perf_counter()
            if verbose:
                print(f'Time taken was {toc - tic}s')

        elif pool.__class__.__name__ is 'MultiPool':
            list(tqdm(
                pool.imap_unordered(rmsythoncut_i, inputs),
                total=count,
                desc='Running RM synth on Stokes I',
                disable=(not verbose)
            )
            )
    else:
        inputs = [[i, args, outdir, freqfile, verbose] for i in range(count)]
        if (pool.__class__.__name__ is 'MPIPool' or
                pool.__class__.__name__ is 'SerialPool'):
            if verbose:
                print('Running RM synth...')
            tic = time.perf_counter()
            list(pool.map(rmsythoncut, inputs))
            toc = time.perf_counter()
            if verbose:
                print(f'Time taken was {toc - tic}s')

        elif pool.__class__.__name__ is 'MultiPool':
            list(tqdm(
                pool.imap_unordered(rmsythoncut, inputs),
                total=count,
                desc='Running RM synth',
                disable=(not verbose)
            )
            )

    if verbose:
        print('Done!')


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
    SPICE-RACS Stage 5:
    Run RMsynthesis on cubelets.

    Note: Runs on brightest sources first.

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

    parser.add_argument("--validate", dest="validate", action="store_true",
                        help="Run on Stokes I [False].")

    parser.add_argument("--limit", dest="limit", default=None,
                        type=int, help="Limit number of sources [All].")

    # RM-tools args
    parser.add_argument("-w", dest="weightType", default="uniform",
                        help="weighting [uniform] (all 1s) or 'variance'.")
    parser.add_argument("-t", dest="fitRMSF", action="store_true",
                        help="Fit a Gaussian to the RMSF [False]")
    parser.add_argument("-l", dest="phiMax_radm2", type=float, default=None,
                        help="Absolute max Faraday depth sampled (overrides NSAMPLES) [Auto].")
    parser.add_argument("-d", dest="dPhi_radm2", type=float, default=None,
                        help="Width of Faraday depth channel [Auto].")
    parser.add_argument("-o", dest="prefixOut", default="",
                        help="Prefix to prepend to output files [None].")
    parser.add_argument("-s", dest="nSamples", type=float, default=5,
                        help="Number of samples across the FWHM RMSF.")
    parser.add_argument("-f", dest="write_seperate_FDF", action="store_false",
                        help="Store different Stokes as FITS extensions [False, store as seperate files].")
    parser.add_argument("-v", dest="verbose", action="store_true",
                        help="Verbose [False].")
    parser.add_argument("-R", dest="not_RMSF", action="store_true",
                        help="Skip calculation of RMSF? [False]")
    parser.add_argument("-rmv", dest="rm_verbose", action="store_true",
                        help="Verbose RMsynth [False].")
    parser.add_argument("-s", dest="NSAMPLES", default=None,
                        type=int, help="Number of samples across the FWHM RMSF.")

    group = parser.add_mutually_exclusive_group()

    group.add_argument("--ncores", dest="n_cores", default=1,
                       type=int, help="Number of processes (uses multiprocessing).")
    group.add_argument("--mpi", dest="mpi", default=False,
                       action="store_true", help="Run with MPI.")

    args = parser.parse_args()

    if args.validate:
        pool = schwimmbad.SerialPool
    else:
        pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)

    verbose = args.verbose

    if args.mpi:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)

    if verbose:
        print(f"Using pool: {pool.__class__.__name__}")

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

    main(pool, args, verbose=verbose)


if __name__ == "__main__":
    cli()
