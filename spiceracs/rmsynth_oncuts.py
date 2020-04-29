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
import pdb
print = functools.partial(print, f'[{psutil.Process().cpu_num()}]', flush=True)


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


def rms_1d(data):
    """Compute RMS from bounding pixels
    """
    Nfreq, Ndec, Nra = data.shape
    mask = np.ones((Ndec, Nra), dtype=np.bool)
    mask[3:-3, 3:-3] = False
    rms = np.nanstd(data[:, mask], axis=1)
    return rms



class ForkedPdb(pdb.Pdb):
    """A Pdb subclass that may be used
    from a forked multiprocessing child

    """
    def interaction(self, *args, **kwargs):
        _stdin = sys.stdin
        try:
            sys.stdin = open('/dev/stdin')
            pdb.Pdb.interaction(self, *args, **kwargs)
        finally:
            sys.stdin = _stdin

def estimate_noise_annulus(x_center, y_center, cube):
    """
    Noise estimation for annulus taken around point source. Annulus has fixed 
    inner radius of 10 and outer radius of 31. Function makes an annulus shaped
    mask, then for each source applies the mask at each frequency and takes the
    standard deviation.

    ​Inputs: Array of sets of pixel coordinates (y-position,x-position) for 
    sources, Stokes cube (assumes 4 axes), array of flagged channels (can be an
    empty array), number of frequency channels.

    ​Output: 2D array of standard deviation values with shape (length of 
    coordinate array, number of unflagged frequency channels).
    """

    inner_radius = 10
    outer_radius = 31

    lenfreq = cube.shape[0]
    naxis = len(cube.shape)
    err = np.zeros(lenfreq)
    try:
        y, x = np.ogrid[-1*outer_radius:outer_radius +
                        1, -1*outer_radius:outer_radius+1]
        grid_mask = np.logical_or(
            x**2+y**2 < inner_radius**2, x**2+y**2 > outer_radius**2)
        for i in range(lenfreq):
            if naxis == 4:
                grid = cube[i, 0, y_center-outer_radius:y_center+outer_radius+1,
                            x_center-outer_radius:x_center+outer_radius+1]
            else:  # naxis ==3
                grid = cube[i, y_center-outer_radius:y_center+outer_radius+1,
                            x_center-outer_radius:x_center+outer_radius+1]

            # Calculate the MADFM, and convert to standard sigma:
            noisepix = np.ma.masked_array(grid, grid_mask)
            err[i] = np.ma.median(np.ma.fabs(
                noisepix - np.ma.median(noisepix))) / 0.6745
            return err
    except np.ma.core.MaskError:
        print('Too many NaNs in cutout - skipping...')
        err *= np.nan
        return err
    

########################################################################
#def estimate_noise_annulus_np(x_center, y_center, cube):
#    inner_radius = 10
#    outer_radius = 31
#    z = np.arange(cube.shape[0])
#    y = np.arange(cube.shape[1])
#    x = np.arange(cube.shape[2])
#    Z, Y, X = np.meshgrid(z, y, x, indexing='ij')
#    R = np.hypot(X-x_center, Y-y_center)
#    mask = np.logical_or(R < inner_radius, R > outer_radius)
#    noise_pix = np.ma.masked_array(cube, mask)
#    err = np.ma.median(
#        np.ma.fabs(
#            noise_pix - np.ma.median(
#                noise_pix, axis=(1, 2)
#            )[:, np.newaxis, np.newaxis]
#        ), axis=(1, 2)
#    ) / 0.6745
#    return err
#######################################################################

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
        #rmsi = rms_1d(dataI)
    rmsi = estimate_noise_annulus(
        dataI.shape[2]//2,
        dataI.shape[1]//2,
        dataI
    )
    rmsi[rmsi == 0] = np.nan
    rmsi[np.isnan(rmsi)] = np.nanmedian(rmsi)

    #rmsq = rms_1d(dataQ)
    rmsq = estimate_noise_annulus(
        dataQ.shape[2]//2,
        dataQ.shape[1]//2,
        dataQ
    )
    rmsq[rmsq == 0] = np.nan
    rmsq[np.isnan(rmsq)] = np.nanmedian(rmsq)

    #rmsu = rms_1d(dataU)
    rmsu = estimate_noise_annulus(
        dataU.shape[2]//2,
        dataU.shape[1]//2,
        dataU
    )
    rmsu[rmsu == 0] = np.nan
    rmsu[np.isnan(rmsu)] = np.nanmedian(rmsu)

    if np.isnan(rmsi).all() or np.isnan(rmsq).all() or np.isnan(rmsu).all():
        return
    if verbose:
        return
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
        import ipdb; ipdb.set_trace()
        qarr = dataQ[:, y, x]
        uarr = dataU[:, y, x]
        iarr = dataI[:, y, x]

        iarr[iarr == 0] = np.nan
        qarr[qarr == 0] = np.nan
        uarr[uarr == 0] = np.nan

        if np.isnan(qarr).all() or np.isnan(uarr).all():
            return
        else:
            if clargs.noStokesI:
                idx = np.isnan(qarr) | np.isnan(uarr)
                data = [np.array(freq), qarr,
                        uarr, rmsq, rmsu]
            else:
                if np.isnan(iarr).all():
                    return
                else:
                    idx = np.isnan(qarr) | np.isnan(uarr) | np.isnan(iarr)
                    data = [np.array(freq), iarr, qarr,
                            uarr, rmsi, rmsq, rmsu]
            # Run 1D RM-synthesis on the spectra

            # if clargs.database:
            #    myquery = {"island_name": iname}
            #    data_f = [np.array(freq), iarr, qarr, uarr, rmsi, rmsq, rmsu]
            #    json_data = json.loads(json.dumps(
            #        {f"comp_{comp+1}_spectra": data}, cls=MyEncoder))
            #    newvalues = {"$set": json_data}
            #    mycol.update_one(myquery, newvalues)
            #    mycol.update_one(myquery, newvalues)

            np.savetxt(f"{prefix}.dat", np.vstack(data).T, delimiter=' ')

            try:
                mDict, aDict = do_RMsynth_1D.run_rmsynth(data=data,
                                                         polyOrd=clargs.polyOrd,
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
                if clargs.savePlots:
                    if verbose:
                        print("Plotting the input data and spectral index fit.")
                    from RMutils.util_plotTk import plot_Ipqu_spectra_fig
                    from RMutils.util_misc import poly5

                    if clargs.noStokesI:
                        IArr = np.ones_like(qarr[~idx])
                        Ierr = np.zeros_like(qarr[~idx])
                    else:
                        IArr = iarr[~idx]
                        Ierr = rmsi[~idx]

                    IModArr, qArr, uArr, dqArr, duArr, fitDict = \
                        create_frac_spectra(freqArr=np.array(freq)[~idx]/1e9,
                                            IArr=IArr,
                                            QArr=qarr[~idx],
                                            UArr=uarr[~idx],
                                            dIArr=Ierr,
                                            dQArr=rmsq[~idx],
                                            dUArr=rmsu[~idx],
                                            polyOrd=clargs.polyOrd,
                                            verbose=False,
                                            debug=False)

                    freqHirArr_Hz = np.linspace(
                        mDict['min_freq'], mDict['max_freq'], 10000)
                    coef = np.array(
                        mDict["polyCoeffs"].split(',')).astype(float)
                    IModHirArr = poly5(coef)(freqHirArr_Hz/1e9)
                    plot_Ipqu_spectra_fig(freqArr_Hz=np.array(freq)[~idx],
                                          IArr=iarr[~idx],
                                          qArr=qArr,
                                          uArr=uArr,
                                          dIArr=Ierr,
                                          dqArr=dqArr,
                                          duArr=duArr,
                                          freqHirArr_Hz=freqHirArr_Hz,
                                          IModArr=IModHirArr,
                                          fig=None,
                                          units='Jy/beam')
                    plotname = f'{outdir}/plots/{iname}_specfig.png'
                    plt.savefig(plotname, dpi=75, bbox_inches='tight')
                do_RMsynth_1D.saveOutput(
                    mDict, aDict, prefix, clargs.rm_verbose)
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

    if args.savePlots:
        plotdir = f'{outdir}/plots'
        try:
            os.mkdir(plotdir)
            print('Made plot directory.')
        except FileExistsError:
            print('Directory exists.')

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
    parser.add_argument("-sp", dest="savePlots", action="store_true",
                        help="save the plots [False].")
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
    parser.add_argument("-o", dest="polyOrd", type=int, default=3,
                        help="polynomial order to fit to I spectrum [3].")
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
