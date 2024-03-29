#!/usr/bin/env python3
"""Run RM-synthesis on cutouts in parallel"""
import os
import sys
import warnings
from glob import glob
from pathlib import Path
from pprint import pformat
from shutil import copyfile
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pymongo
from prefect import flow, task, unmapped
from RMtools_1D import do_RMclean_1D
from RMtools_3D import do_RMclean_3D

from arrakis.logger import logger
from arrakis.utils.database import get_db, test_db
from arrakis.utils.pipeline import logo_str


@task(name="1D RM-CLEAN")
def rmclean1d(
    comp: dict,
    outdir: str,
    cutoff: float = -3,
    maxIter=10000,
    gain=0.1,
    showPlots=False,
    savePlots=False,
    rm_verbose=True,
    window=None,
) -> pymongo.UpdateOne:
    """1D RM-CLEAN

    Args:
        comp (dict): Mongo entry for component.
        outdir (str): Output directory.
        cutoff (float, optional): CLEAN cutouff (in sigma). Defaults to -3.
        maxIter (int, optional): Maximum CLEAN interation. Defaults to 10000.
        gain (float, optional): CLEAN gain. Defaults to 0.1.
        showPlots (bool, optional): Show CLEAN plots. Defaults to False.
        savePlots (bool, optional): Save CLEAN plots. Defaults to False.
        rm_verbose (bool, optional): Verbose RM-CLEAN. Defaults to True.

    Returns:
        pymongo.UpdateOne: MongoDB update query.
    """
    iname = comp["Source_ID"]
    cname = comp["Gaussian_ID"]

    logger.debug(f"Working on {comp}")
    try:
        rm1dfiles = comp["rm1dfiles"]
        fdfFile = os.path.join(outdir, f"{rm1dfiles['FDF_dirty']}")
        rmsfFile = os.path.join(outdir, f"{rm1dfiles['RMSF']}")
        weightFile = os.path.join(outdir, f"{rm1dfiles['weights']}")
        rmSynthFile = os.path.join(outdir, f"{rm1dfiles['summary_json']}")

        prefix = os.path.join(os.path.abspath(os.path.dirname(fdfFile)), cname)

        # Sanity checks
        for f in [weightFile, fdfFile, rmsfFile, rmSynthFile]:
            logger.debug(f"Checking {os.path.abspath(f)}")
            if not os.path.exists(f):
                logger.fatal("File does not exist: '{:}'.".format(f))
                sys.exit()
        nBits = 32
        mDict, aDict = do_RMclean_1D.readFiles(
            fdfFile, rmsfFile, weightFile, rmSynthFile, nBits
        )

        # Run RM-CLEAN on the spectrum
        outdict, arrdict = do_RMclean_1D.run_rmclean(
            mDict=mDict,
            aDict=aDict,
            cutoff=cutoff,
            maxIter=maxIter,
            gain=gain,
            nBits=nBits,
            showPlots=showPlots,
            verbose=rm_verbose,
            prefixOut=prefix,
            saveFigures=savePlots,
            window=window,
        )
        # Ensure JSON serializable
        for k, v in outdict.items():
            if isinstance(v, np.float_):
                outdict[k] = float(v)
            elif isinstance(v, np.float32):
                outdict[k] = float(v)
            elif isinstance(v, np.int_):
                outdict[k] = int(v)
            elif isinstance(v, np.int32):
                outdict[k] = int(v)
            elif isinstance(v, np.ndarray):
                outdict[k] = v.tolist()

        # Save output
        do_RMclean_1D.saveOutput(outdict, arrdict, prefixOut=prefix, verbose=rm_verbose)
        if savePlots:
            plt.close("all")
            plotdir = os.path.join(outdir, "plots")
            plot_files = glob(
                os.path.join(os.path.abspath(os.path.dirname(fdfFile)), "*.pdf")
            )
            for src in plot_files:
                base = os.path.basename(src)
                dst = os.path.join(plotdir, base)
                copyfile(src, dst)
        # Load into Mongo
        myquery = {"Gaussian_ID": cname}

        newvalues = {
            "$set": {
                "rmclean1d": True,
                "rmclean_summary": outdict,
            },
        }
    except KeyError:
        logger.critical("Failed to load data! RM-CLEAN not applied to component!")
        logger.critical(f"Island is {iname}, component is {cname}")
        myquery = {"Gaussian_ID": cname}

        newvalues = {
            "$set": {
                "rmclean1d": False,
            },
        }
    return pymongo.UpdateOne(myquery, newvalues)


@task(name="3D RM-CLEAN")
def rmclean3d(
    island: dict,
    outdir: str,
    cutoff: float = -3,
    maxIter=10000,
    gain=0.1,
    rm_verbose=False,
) -> pymongo.UpdateOne:
    """Run RM-CLEAN on 3D cube

    Args:
        island (dict): MongoDB island entry.
        outdir (str): Output directory.
        cutoff (float, optional): CLEAN cutoff (in sigma). Defaults to -3.
        maxIter (int, optional): Max CLEAN iterations. Defaults to 10000.
        gain (float, optional): CLEAN gain. Defaults to 0.1.
        rm_verbose (bool, optional): Verbose output. Defaults to False.

    Returns:
        pymongo.UpdateOne: MongoDB update query.
    """
    """3D RM-CLEAN

    Args:
        island_id (str): RACS Island ID
        host (str): MongoDB host
        field (str): RACS field
        cutoff (int, optional): CLEAN cutoff. Defaults to -3.
        maxIter (int, optional): CLEAN max iterations. Defaults to 10000.
        gain (float, optional): CLEAN gain. Defaults to 0.1.
        rm_verbose (bool, optional): Verbose RM-CLEAN. Defaults to False.
    """
    iname = island["Source_ID"]
    prefix = f"{iname}_"
    rm3dfiles = island["rm3dfiles"]

    cleanFDF, ccArr, iterCountArr, residFDF, headtemp = do_RMclean_3D.run_rmclean(
        fitsFDF=os.path.join(outdir, rm3dfiles["FDF_real_dirty"]),
        fitsRMSF=os.path.join(outdir, rm3dfiles["RMSF_tot"]),
        cutoff=cutoff,
        maxIter=maxIter,
        gain=gain,
        chunksize=None,
        nBits=32,
        verbose=rm_verbose,
    )

    # Write results to disk
    do_RMclean_3D.writefits(
        cleanFDF,
        ccArr,
        iterCountArr,
        residFDF,
        headtemp,
        prefixOut=prefix,
        outDir=os.path.abspath(
            os.path.dirname(os.path.join(outdir, rm3dfiles["FDF_real_dirty"]))
        ),
        write_separate_FDF=True,
        verbose=rm_verbose,
    )
    # Load into Mongo
    myquery = {"Source_ID": iname}
    newvalues = {"$set": {"rmclean3d": True}}
    return pymongo.UpdateOne(myquery, newvalues)


@flow(name="RM-CLEAN on cutouts")
def main(
    field: str,
    outdir: Path,
    host: str,
    epoch: int,
    username: Optional[str] = None,
    password: Optional[str] = None,
    dimension="1d",
    database=False,
    savePlots=True,
    limit: Optional[int] = None,
    cutoff: float = -3,
    maxIter=10000,
    gain=0.1,
    window=None,
    showPlots=False,
    rm_verbose=False,
):
    """Main script

    Args:
        field (str): RACS field name.
        outdir (str): Output directory.
        host (str): MongoDB host IP.
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo password. Defaults to None.
        dimension (str, optional): Which dimension to run RM-CLEAN. Defaults to "1d".
        verbose (bool, optional): Verbose output. Defaults to True.
        database (bool, optional): Update database. Defaults to False.
        savePlots (bool, optional): Save plots. Defaults to True.
        validate (bool, optional): Run validation. Defaults to False.
        limit (int, optional): Limit number of sources processed. Defaults to None.
        cutoff (float, optional): CLEAN cutoff (in sigma). Defaults to -3.
        maxIter (int, optional): Max CLEAN iterations. Defaults to 10000.
        gain (float, optional): Clean gain. Defaults to 0.1.
        showPlots (bool, optional): Show interactive plots. Defaults to False.
        rm_verbose (bool, optional): Verbose output from RM-CLEAN. Defaults to False.
    """
    outdir = os.path.abspath(outdir)
    outdir = os.path.join(outdir, "cutouts")

    # default connection (ie, local)
    beams_col, island_col, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )

    query = {"$and": [{f"beams.{field}": {"$exists": True}}]}

    beams = list(beams_col.find(query).sort("Source_ID"))
    all_island_ids = sorted(beams_col.distinct("Source_ID", query))

    if dimension == "3d":
        query = {"$and": [{"Source_ID": {"$in": all_island_ids}}, {"rmsynth3d": True}]}

        islands = list(
            island_col.find(
                query,
                # Only get required values
                {
                    "Source_ID": 1,
                    "rm3dfiles": 1,
                },
            ).sort("Source_ID")
        )
        island_ids = [doc["Source_ID"] for doc in islands]
        n_island = island_col.count_documents(query)
        island_col.update(query, {"$set": {"rmclean3d": False}})

    elif dimension == "1d":
        query = {"$and": [{"Source_ID": {"$in": all_island_ids}}, {"rmsynth1d": True}]}

        components = list(
            comp_col.find(
                query,
                # Only get required values
                {
                    "Source_ID": 1,
                    "Gaussian_ID": 1,
                    "rm1dfiles": 1,
                },
            ).sort("Source_ID")
        )
        n_comp = comp_col.count_documents(query)
        comp_col.update_many(query, {"$set": {"rmclean1d": False}})

    if limit is not None:
        count = limit
        n_comp = count
        n_island = count
        # component_ids = component_ids[:count]

    if dimension == "1d":
        logger.info(f"Running RM-CLEAN on {n_comp} components")
        outputs = rmclean1d.map(
            comp=components,
            outdir=unmapped(outdir),
            cutoff=unmapped(cutoff),
            maxIter=unmapped(maxIter),
            gain=unmapped(gain),
            showPlots=unmapped(showPlots),
            savePlots=unmapped(savePlots),
            rm_verbose=unmapped(rm_verbose),
            window=unmapped(window),
        )
    elif dimension == "3d":
        logger.info(f"Running RM-CLEAN on {n_island} islands")

        outputs = rmclean3d.map(
            island=islands,
            outdir=unmapped(outdir),
            cutoff=unmapped(cutoff),
            maxIter=unmapped(maxIter),
            gain=unmapped(gain),
            rm_verbose=unmapped(rm_verbose),
        )

    else:
        raise ValueError(f"Dimension {dimension} not supported.")

    if database:
        logger.info("Updating database...")
        updates = [f.result() for f in outputs]
        if dimension == "1d":
            db_res = comp_col.bulk_write(updates, ordered=False)
            logger.info(pformat(db_res.bulk_api_result))
        elif dimension == "3d":
            db_res = island_col.bulk_write(updates, ordered=False)
            logger.info(pformat(db_res.bulk_api_result))
    logger.info("RM-CLEAN done!")


def cli():
    """Command-line interface"""
    import argparse

    from astropy.utils.exceptions import AstropyWarning

    warnings.simplefilter("ignore", category=AstropyWarning)
    from astropy.io.fits.verify import VerifyWarning

    warnings.simplefilter("ignore", category=VerifyWarning)
    # Help string to be shown using the -h option
    descStr = f"""
    {logo_str}
    Arrakis Stage 6:
    Run RM-CLEAN on cubelets.

    Note: Runs on brightest sources first.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field to mosaic - e.g. 2132-50A."
    )
    parser.add_argument(
        "outdir",
        metavar="outdir",
        type=Path,
        help="Directory containing cutouts (in subdir outdir/cutouts).",
    )

    parser.add_argument(
        "host",
        metavar="host",
        type=str,
        help="Host of mongodb (probably $hostname -i).",
    )

    parser.add_argument(
        "-e",
        "--epoch",
        type=int,
        default=0,
        help="Epoch of observation.",
    )

    parser.add_argument(
        "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )

    parser.add_argument(
        "--dimension",
        dest="dimension",
        default="1d",
        help="How many dimensions for RMsynth [1d] or '3d'.",
    )

    parser.add_argument(
        "-v", dest="verbose", action="store_true", help="verbose output [False]."
    )

    parser.add_argument(
        "-m", dest="database", action="store_true", help="Add data to MongoDB [False]."
    )
    parser.add_argument(
        "-sp", "--savePlots", action="store_true", help="save the plots [False]."
    )

    parser.add_argument(
        "--limit",
        dest="limit",
        default=None,
        type=int,
        help="Limit number of sources [All].",
    )

    # RM-tools args
    parser.add_argument(
        "-c",
        dest="cutoff",
        type=float,
        default=-3,
        help="CLEAN cutoff (+ve = absolute, -ve = sigma) [-3].",
    )
    parser.add_argument(
        "-n",
        dest="maxIter",
        type=int,
        default=10000,
        help="maximum number of CLEAN iterations [10000].",
    )
    parser.add_argument(
        "-g", dest="gain", type=float, default=0.1, help="CLEAN loop gain [0.1]."
    )
    parser.add_argument(
        "-w",
        dest="window",
        type=float,
        default=None,
        help="Further CLEAN in mask to this threshold [False].",
    )
    parser.add_argument(
        "-p", dest="showPlots", action="store_true", help="show the plots [False]."
    )
    parser.add_argument(
        "-rmv", dest="rm_verbose", action="store_true", help="Verbose RM-CLEAN [False]."
    )

    args = parser.parse_args()

    verbose = args.verbose
    rmv = args.rm_verbose
    host = args.host
    test_db(
        host=args.host, username=args.username, password=args.password, verbose=verbose
    )

    if rmv:
        logger.setLevel(
            level=logger.DEBUG,
        )
    elif verbose:
        logger.setLevel(
            level=logger.INFO,
        )
    main(
        field=args.field,
        outdir=Path(args.outdir),
        host=host,
        username=args.username,
        password=args.password,
        dimension=args.dimension,
        database=args.database,
        savePlots=args.savePlots,
        limit=args.limit,
        cutoff=args.cutoff,
        maxIter=args.maxIter,
        gain=args.gain,
        window=args.window,
        showPlots=args.showPlots,
        rm_verbose=args.rm_verbose,
    )


if __name__ == "__main__":
    cli()
