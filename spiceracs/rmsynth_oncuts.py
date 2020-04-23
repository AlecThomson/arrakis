#!/usr/bin/env python
from spiceracs.utils import getfreq, MyEncoder
import json
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
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from RMutils.util_misc import create_frac_spectra
import functools
import psutil
print = functools.partial(print,f'[{psutil.Process().cpu_num()}]', flush=True)


def rmsythoncut3d(args):
    i, clargs, outdir, freq, freqfile, verbose = args

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
        rms = np.nanstd(np.squeeze(hdulist[0].data), axis=(1, 2)) * 3

    prefix = iname

    header, dataQ = do_RMsynth_3D.readFitsCube(qfile, clargs.rm_verbose)

    rmsArr = np.ones_like(dataQ)*rms[:, np.newaxis, np.newaxis]

    # Run 3D RM-synthesis on the cubes
    dataArr = do_RMsynth_3D.run_rmsynth(
        dataQ=np.squeeze(dataQ),
        dataU=np.squeeze(do_RMsynth_3D.readFitsCube(
            ufile, clargs.rm_verbose)[1]),
        freqArr_Hz=do_RMsynth_3D.readFreqFile(freqfile, clargs.rm_verbose),
        dataI=None,
        rmsArr=rmsArr,
        phiMax_radm2=clargs.phiMax_radm2,
        dPhi_radm2=clargs.dPhi_radm2,
        nSamples=clargs.nSamples,
        weightType=clargs.weightType,
        fitRMSF=clargs.fitRMSF,
        nBits=32,
        verbose=clargs.rm_verbose,
        not_rmsf=clargs.not_RMSF)

    # Write to files
    do_RMsynth_3D.writefits(dataArr,
                            headtemplate=header,
                            fitRMSF=clargs.fitRMSF,
                            prefixOut=prefix,
                            outDir=outdir,
                            write_seperate_FDF=True,
                            not_rmsf=clargs.not_RMSF,
                            nBits=32,
                            verbose=clargs.rm_verbose)

    if clargs.database:
        myquery = {"island_name": iname}

        newvalues = {"$set": {"rmsynth3d": True}}
        mycol.update_one(myquery, newvalues)

        newvalues = {"$set": {"rm3dfiles": {
            "FDF_real_dirty": f"{prefix}FDF_real_dirty.fits",
            "FDF_im_dirty": f"{prefix}FDF_im_dirty.fits",
            "FDF_tot_dirty": f"{prefix}FDF_tot_dirty.fits",
            "RMSF_real": f"{prefix}RMSF_real.fits",
            "RMSF_tot": f"{prefix}RMSF_tot.fits",
            "RMSF_FWHM": f"{prefix}RMSF_FWHM.fits"
        }}}
        mycol.update_one(myquery, newvalues)


def rmsythoncut1d(args):
    i, clargs, outdir, freq, freqfile, verbose = args

    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

    # Basic querey
    if clargs.pol and not clargs.unres:
        myquery = {"polarized": True}
    elif clargs.unres and not clargs.pol:
        myquery = {"resolved": False}
    elif clargs.pol and clargs.unres:
        myquery = {"$and": [{"resolved": False}, {"polarized": True}]}

    elif clargs.pol and not clargs.loners:
        myquery = {"polarized": True}
    elif clargs.loners and not clargs.pol:
        myquery = {"n_components": 1}
    elif clargs.pol and clargs.loners:
        myquery = {"$and": [{"n_components": 1}, {"polarized": True}]}
    else:
        myquery = {}

    doc = mycol.find(myquery).sort("flux_peak", -1)
    iname = doc[i]['island_name']
    ifile = f"{outdir}/{doc[i]['i_file']}"
    qfile = f"{outdir}/{doc[i]['q_file']}"
    ufile = f"{outdir}/{doc[i]['u_file']}"
    vfile = f"{outdir}/{doc[i]['v_file']}"

    # with fits.open(vfile) as hdulist:
    #    rms = np.nanstd(np.squeeze(hdulist[0].data), axis=(1, 2)) * 3
    #rms[rms == 0] = np.nan
    # get edge values

    header, dataQ = do_RMsynth_3D.readFitsCube(qfile, clargs.rm_verbose)
    header, dataU = do_RMsynth_3D.readFitsCube(ufile, clargs.rm_verbose)
    header, dataI = do_RMsynth_3D.readFitsCube(ifile, clargs.rm_verbose)

    dataQ = np.squeeze(dataQ)
    dataU = np.squeeze(dataU)
    dataI = np.squeeze(dataI)

    if np.isnan(dataI).all() or np.isnan(dataQ).all() or np.isnan(dataU).all():
        return

    vari = np.nansum([np.nanvar(dataI[:, :, 0:3], axis=(1, 2)),
                      np.nanvar(dataI[:, :, -4:-1], axis=(1, 2)), 
                      np.nanvar(dataI[:, 0:3, :], axis=(1, 2)),
                      np.nanvar(dataI[:, -4:-1, :], axis=(1, 2))], 
                      axis=0)
    vari[vari == 0] = np.nan
    rmsi = np.sqrt(vari)
    rmsi[np.isnan(rmsi)] = np.nanmedian(rmsi)

    varq = np.nansum([np.nanvar(dataQ[:, :, 0:3], axis=(1, 2)),
                      np.nanvar(dataQ[:, :, -4:-1], axis=(1, 2)), 
                      np.nanvar(dataQ[:, 0:3, :], axis=(1, 2)),
                      np.nanvar(dataQ[:, -4:-1, :], axis=(1, 2))], 
                      axis=0)
    varq[varq == 0] = np.nan
    rmsq = np.sqrt(varq)
    rmsq[np.isnan(rmsq)] = np.nanmedian(rmsq)

    varu = np.nansum([np.nanvar(dataU[:, :, 0:3], axis=(1, 2)),
                      np.nanvar(dataU[:, :, -4:-1], axis=(1, 2)), 
                      np.nanvar(dataU[:, 0:3, :], axis=(1, 2)),
                      np.nanvar(dataU[:, -4:-1, :], axis=(1, 2))], 
                      axis=0)
    varu[varu == 0] = np.nan
    rmsu = np.sqrt(varu)
    rmsu[np.isnan(rmsu)] = np.nanmedian(rmsu)

    for comp in range(int(doc[i]['n_components'])):
        if clargs.rm_verbose:
            print(f'Working on component {comp+1}')
        cname = doc[i][f'component_{comp+1}']['component_name']
        prefix = f'{outdir}/{cname}'

        ra = doc[i][f'component_{comp+1}']['ra_deg_cont']
        dec = doc[i][f'component_{comp+1}']['dec_deg_cont']
        wcs = WCS(doc[i]['header']).dropaxis(2)

        x, y, z = np.array(wcs.all_world2pix(
            ra, dec, np.nanmean(freq), 0)).round().astype(int)

        qarr = np.nansum(dataQ[:, y-1:y+1+1, x-1:x+1+1], axis=(1, 2))
        uarr = np.nansum(dataU[:, y-1:y+1+1, x-1:x+1+1], axis=(1, 2))
        iarr = np.nansum(dataI[:, y-1:y+1+1, x-1:x+1+1], axis=(1, 2))

        iarr[iarr == 0] = np.nan
        qarr[qarr == 0] = np.nan
        uarr[uarr == 0] = np.nan

        if np.isnan(qarr).all() or np.isnan(uarr).all():
            return
        else:
            if clargs.noStokesI:
                idx = np.isnan(qarr) | np.isnan(uarr)
                data = [np.array(freq)[~idx], qarr[~idx], uarr[~idx], rmsq[~idx], rmsu[~idx]]
            else:
                if np.isnan(iarr).all():
                    return
                else:
                    idx = np.isnan(qarr) | np.isnan(uarr) | np.isnan(iarr)
                    data = [np.array(freq)[~idx], iarr[~idx], qarr[~idx], uarr[~idx], rmsi[~idx], rmsq[~idx], rmsu[~idx]]
            # Run 1D RM-synthesis on the spectra

            np.savetxt(f"{prefix}.dat", np.vstack(data).T, delimiter=' ')

            try:
                mDict, aDict = do_RMsynth_1D.run_rmsynth(data=data,
                                                        phiMax_radm2=clargs.phiMax_radm2,
                                                        dPhi_radm2=clargs.dPhi_radm2,
                                                        nSamples=clargs.nSamples,
                                                        weightType=clargs.weightType,
                                                        fitRMSF=clargs.fitRMSF,
                                                        noStokesI=clargs.noStokesI,
                                                        nBits=32,
                                                        showPlots=clargs.showPlots,
                                                        verbose=clargs.rm_verbose,
                                                        debug=clargs.debug)
                import ipdb; ipdb.set_trace()
                if clargs.savePlots:
                    if verbose: log("Plotting the input data and spectral index fit.")
                    from RMutils.util_plotTk import plot_Ipqu_spectra_fig
                    from RMutils.util_misc import poly5
                    freqHirArr_Hz =  np.linspace(freqArr_Hz[0], freqArr_Hz[-1], 10000)
                    coef = np.array(mDict["polyCoeffs"].split(',')).astype(float)
                    IModHirArr = poly5(coef)(freqHirArr_Hz/1e9)
                    specFig = plt.figure(figsize=(12.0, 8))
                    plot_Ipqu_spectra_fig(freqArr_Hz     = freqArr_Hz,
                                        IArr           = IArr,
                                        qArr           = qArr,
                                        uArr           = uArr,
                                        dIArr          = dIArr,
                                        dqArr          = dqArr,
                                        duArr          = duArr,
                                        freqHirArr_Hz  = freqHirArr_Hz,
                                        IModArr        = IModHirArr,
                                        fig            = specFig,
                                        units          = units)

                do_RMsynth_1D.saveOutput(mDict, aDict, prefix, clargs.rm_verbose)
            except ValueError:
                return
            if clargs.database:
                myquery = {"island_name": iname}

                newvalues = {"$set": {f"rm1dfiles_comp_{comp+1}": {
                    "FDF_dirty": f"{cname}_FDFdirty.dat",
                    "RMSF": f"{cname}_RMSF.dat",
                    "weights": f"{cname}_weight.dat",
                    "summary_dat": f"{cname}_RMsynth.dat",
                    "summary_json": f"{cname}_RMsynth.json",
                }}}
                mycol.update_one(myquery, newvalues)

                myquery = {"island_name": iname}
                newvalues = {"$set": {f"comp_{comp+1}_rmsynth1d": True}}
                mycol.update_one(myquery, newvalues)

    if clargs.database:
        # Load into Mongo
        myquery = {"island_name": iname}

        newvalues = {"$set": {f"rmsynth1d": True}}
        mycol.update_one(myquery, newvalues)


def rmsythoncut_i(args):
    i, clargs, freq, outdir, verbose = args

    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

    # Basic querey
    if clargs.pol and not clargs.unres:
        myquery = {"polarized": True}
    elif clargs.unres and not clargs.pol:
        myquery = {"resolved": False}
    elif clargs.pol and clargs.unres:
        myquery = {"$and": [{"resolved": False}, {"polarized": True}]}

    elif clargs.pol and not clargs.loners:
        myquery = {"polarized": True}
    elif clargs.loners and not clargs.pol:
        myquery = {"n_components": 1}
    elif clargs.pol and clargs.loners:
        myquery = {"$and": [{"n_components": 1}, {"polarized": True}]}
    else:
        myquery = {}

    doc = mycol.find(myquery).sort("flux_peak", -1)

    iname = doc[i]['island_name']
    data = np.array(SpectralCube.read(f"{outdir}/{doc[i]['i_file']}"))

    for comp in range(int(doc[i]['n_components'])):
        cname = doc[i][f'component_{comp+1}']['component_name']
        prefix = f'{outdir}/{cname}'
        # Get source peak from Selavy
        ra = doc[i]['ra_deg_cont']
        dec = doc[i]['dec_deg_cont']
        wcs = WCS(doc[i]['header'])

        x, y, z = np.array(wcs.all_world2pix(
            ra, dec, np.nanmean(freq), 0)).round().astype(int)

        mom = np.nansum(data, axis=0)

        plt.ion()
        plt.figure()
        plt.imshow(mom, origin='lower', cmap='cubehelix_r')
        plt.scatter(x, y, c='r', marker='x')
        plt.show()
        _ = input("Press [enter] to continue")  # wait for input from the user
        plt.close()    # close the figure to show the next one.

        data = np.nansum(data[:, y-1:y+1+1, x-1:x+1+1], axis=(1, 2))

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
        do_RMsynth_1D.saveOutput(mDict, aDict, prefix, verbose=verbose)


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
    elif args.unres and not args.pol:
        myquery = {"resolved": False}
    elif args.pol and args.unres:
        myquery = {"$and": [{"resolved": False}, {"polarized": True}]}

    elif args.pol and not args.loners:
        myquery = {"polarized": True}
    elif args.loners and not args.pol:
        myquery = {"n_components": 1}
    elif args.pol and args.loners:
        myquery = {"$and": [{"n_components": 1}, {"polarized": True}]}

    else:
        myquery = {}

    mydoc = mycol.find(myquery).sort("flux_peak", -1)
    count = mycol.count_documents(myquery)

    if args.limit is not None:
        count = args.limit

    # Make frequency file
    freq, freqfile = getfreq(
        f"{outdir}/{mydoc[0]['q_file']}", outdir=outdir, filename='frequencies.txt', verbose=verbose)
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
    elif args.dimension == '1d':
        inputs = [[i, args, outdir, freq, freqfile, verbose]
                  for i in range(count)]
        if (pool.__class__.__name__ is 'MPIPool' or
                pool.__class__.__name__ is 'SerialPool'):
            if verbose:
                print('Running 1D RM synth...')
            tic = time.perf_counter()
            list(pool.map(rmsythoncut1d, inputs))
            toc = time.perf_counter()
            if verbose:
                print(f'Time taken was {toc - tic}s')

        elif pool.__class__.__name__ is 'MultiPool':
            list(tqdm(
                pool.imap_unordered(rmsythoncut1d, inputs),
                total=count,
                desc='Running 1D RM synth',
                disable=(not verbose)
            )
            )

    elif args.dimension == '3d':
        inputs = [[i, args, outdir, freq, freqfile, verbose]
                  for i in range(count)]
        if (pool.__class__.__name__ is 'MPIPool' or
                pool.__class__.__name__ is 'SerialPool'):
            if verbose:
                print('Running 3D RM synth...')
            tic = time.perf_counter()
            list(pool.map(rmsythoncut3d, inputs))
            toc = time.perf_counter()
            if verbose:
                print(f'Time taken was {toc - tic}s')

        elif pool.__class__.__name__ is 'MultiPool':
            list(tqdm(
                pool.imap_unordered(rmsythoncut3d, inputs),
                total=count,
                desc='Running 3D RM synth',
                disable=(not verbose)
            )
            )

    pool.close()

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
    warnings.simplefilter('ignore', category=RuntimeWarning)
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

    parser.add_argument("--dimension", dest="dimension", default="1d",
                        help="How many dimensions for RMsynth [1d] or '3d'.")

    parser.add_argument("-v", dest="verbose", action="store_true",
                        help="verbose output [False].")

    parser.add_argument(
        "-m",
        dest="database",
        action="store_true",
        help="Add data to MongoDB [False]."
    )

    parser.add_argument("--pol", dest="pol", action="store_true",
                        help="Run on polarized sources [False].")

    parser.add_argument("--unres", dest="unres", action="store_true",
                        help="Run on unresolved sources [False].")

    parser.add_argument("--validate", dest="validate", action="store_true",
                        help="Run on Stokes I [False].")

    parser.add_argument("--limit", dest="limit", default=None,
                        type=int, help="Limit number of sources [All].")

    parser.add_argument("--loners", dest="loners", action="store_true",
                        help="Run on single component sources [False].")

    # RM-tools args
    parser.add_argument("-w", dest="weightType", default="uniform",
                        help="weighting [uniform] (all 1s) or 'variance'.")
    parser.add_argument("-t", dest="fitRMSF", action="store_true",
                        help="Fit a Gaussian to the RMSF [False]")
    parser.add_argument("-l", dest="phiMax_radm2", type=float, default=None,
                        help="Absolute max Faraday depth sampled (overrides NSAMPLES) [Auto].")
    parser.add_argument("-d", dest="dPhi_radm2", type=float, default=None,
                        help="Width of Faraday depth channel [Auto].")
    parser.add_argument("-s", dest="nSamples", type=float, default=5,
                        help="Number of samples across the FWHM RMSF.")
    parser.add_argument("-o", dest="polyOrd", type=int, default=2,
                        help="polynomial order to fit to I spectrum [2].")
    parser.add_argument("-i", dest="noStokesI", action="store_true",
                        help="ignore the Stokes I spectrum [False].")
    parser.add_argument("-p", dest="showPlots", action="store_true",
                        help="show the plots [False].")
    parser.add_argument("-R", dest="not_RMSF", action="store_true",
                        help="Skip calculation of RMSF? [False]")
    parser.add_argument("-rmv", dest="rm_verbose", action="store_true",
                        help="Verbose RMsynth [False].")
    parser.add_argument("-D", dest="debug", action="store_true",
                        help="turn on debugging messages & plots [False].")

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

    main(pool, args, verbose=verbose)


if __name__ == "__main__":
    cli()
