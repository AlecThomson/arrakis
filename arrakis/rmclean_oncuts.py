#!/usr/bin/env python3
"""Run RM-synthesis on cutouts in parallel"""
import argparse
import logging
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

from arrakis import rmsynth_oncuts
from arrakis.logger import UltimateHelpFormatter, logger
from arrakis.utils.database import get_db, test_db
from arrakis.utils.pipeline import generic_parser, logo_str


@task(name="1D RM-CLEAN")
def rmclean1d(
    field: str,
    comp: dict,
    outdir: Path,
    cutoff: float = -3,
    maxIter=10000,
    gain=0.1,
    sbid: Optional[int] = None,
    showPlots=False,
    savePlots=False,
    rm_verbose=True,
    window=None,
) -> pymongo.UpdateOne:
    """1D RM-CLEAN

    Args:
        field (str): RACS field name.
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
    save_name = field if sbid is None else f"{field}_{sbid}"
    try:
        rm1dfiles = comp["rm1dfiles"]
        fdfFile = outdir / f"{rm1dfiles['FDF_dirty']}"
        rmsfFile = outdir / f"{rm1dfiles['RMSF']}"
        weightFile = outdir / f"{rm1dfiles['weights']}"
        rmSynthFile = outdir / f"{rm1dfiles['summary_json']}"

        prefix = os.path.join(os.path.abspath(os.path.dirname(fdfFile)), cname)

        # Sanity checks
        for f in [weightFile, fdfFile, rmsfFile, rmSynthFile]:
            logger.debug(f"Checking {f.absolute()}")
            if not f.exists():
                logger.fatal(f"File does not exist: '{f}'.")
                raise FileNotFoundError(f"File does not exist: '{f}'")

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
            plotdir = outdir / "plots"
            plot_files = list(fdfFile.parent.glob("*.pdf"))
            for plot_file in plot_files:
                copyfile(plot_file, plotdir / plot_file.name)

        # Load into Mongo
        myquery = {"Gaussian_ID": cname}

        newvalues = {
            "$set": {
                save_name: {
                    "rmclean1d": True,
                    "rmclean_summary": outdict,
                },
            }
        }
    except KeyError:
        logger.critical("Failed to load data! RM-CLEAN not applied to component!")
        logger.critical(f"Island is {iname}, component is {cname}")
        myquery = {"Gaussian_ID": cname}

        newvalues = {
            "$set": {
                save_name: {
                    "rmclean1d": False,
                },
            }
        }
    return pymongo.UpdateOne(myquery, newvalues)


@task(name="3D RM-CLEAN")
def rmclean3d(
    field: str,
    island: dict,
    outdir: Path,
    sbid: Optional[int] = None,
    cutoff: float = -3,
    maxIter=10000,
    gain=0.1,
    rm_verbose=False,
) -> pymongo.UpdateOne:
    """Run RM-CLEAN on 3D cube

    Args:
        island (dict): MongoDB island entry.
        outdir (Path): Output directory.
        cutoff (float, optional): CLEAN cutoff (in sigma). Defaults to -3.
        maxIter (int, optional): Max CLEAN iterations. Defaults to 10000.
        gain (float, optional): CLEAN gain. Defaults to 0.1.
        rm_verbose (bool, optional): Verbose output. Defaults to False.

    Returns:
        pymongo.UpdateOne: MongoDB update query.
    """

    iname = island["Source_ID"]
    prefix = f"{iname}_"
    rm3dfiles = island["rm3dfiles"]

    cleanFDF, ccArr, iterCountArr, residFDF, headtemp = do_RMclean_3D.run_rmclean(
        fitsFDF=(outdir / rm3dfiles["FDF_real_dirty"]).as_posix(),
        fitsRMSF=(outdir / rm3dfiles["RMSF_tot"]).as_posix(),
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
        outDir=(outdir / rm3dfiles["FDF_real_dirty"]).parent.absolute().as_posix(),
        write_separate_FDF=True,
        verbose=rm_verbose,
    )
    # Load into Mongo
    save_name = field if sbid is None else f"{field}_{sbid}"
    myquery = {"Source_ID": iname}
    newvalues = {"$set": {save_name: {"rmclean3d": True}}}
    return pymongo.UpdateOne(myquery, newvalues)


@flow(name="RM-CLEAN on cutouts")
def main(
    field: str,
    outdir: Path,
    host: str,
    epoch: int,
    sbid: Optional[int] = None,
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
    """Run RM-CLEAN on cutouts flow

    Args:
        field (str): RACS field name.
        outdir (Path): Output directory.
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
    outdir = outdir.absolute() / "cutouts"

    # default connection (ie, local)
    beams_col, island_col, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )

    query = {"$and": [{f"beams.{field}": {"$exists": True}}]}
    if sbid is not None:
        query["$and"].append({f"beams.{field}.SBIDs": sbid})

    all_island_ids = sorted(beams_col.distinct("Source_ID", query))

    if dimension == "3d":
        query = {
            "$and": [
                {"Source_ID": {"$in": all_island_ids}},
                {
                    (
                        f"{field}.rmsynth3d"
                        if sbid is None
                        else f"{field}_{sbid}.rmsynth3d"
                    ): True
                },
            ]
        }

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
        n_island = island_col.count_documents(query)
        island_col.update(
            query,
            {
                "$set": {
                    (
                        f"{field}.rmclean3d"
                        if sbid is None
                        else f"{field}_{sbid}.rmclean3d"
                    ): False
                }
            },
        )

    elif dimension == "1d":
        query = {
            "$and": [
                {"Source_ID": {"$in": all_island_ids}},
                {
                    (
                        f"{field}.rmsynth1d"
                        if sbid is None
                        else f"{field}_{sbid}.rmsynth1d"
                    ): True
                },
            ]
        }

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
        comp_col.update_many(
            query,
            {
                "$set": {
                    (
                        f"{field}.rmclean1d"
                        if sbid is None
                        else f"{field}_{sbid}.rmclean1d"
                    ): True
                }
            },
        )

    if limit is not None:
        count = limit
        n_comp = count
        n_island = count

    if dimension == "1d":
        logger.info(f"Running RM-CLEAN on {n_comp} components")
        outputs = rmclean1d.map(
            comp=components,
            field=unmapped(field),
            sbid=unmapped(sbid),
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
            field=unmapped(field),
            island=islands,
            sbid=unmapped(sbid),
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


def clean_parser(parent_parser: bool = False) -> argparse.ArgumentParser:
    # Help string to be shown using the -h option
    descStr = f"""
    {logo_str}
    Arrakis Stage 6:
    Run RM-CLEAN on cubelets.

    Note: Runs on brightest sources first.

    """

    # Parse the command line options
    clean_parser = argparse.ArgumentParser(
        add_help=not parent_parser,
        description=descStr,
        formatter_class=UltimateHelpFormatter,
    )
    parser = clean_parser.add_argument_group("rm-clean arguments")

    # RM-tools args
    parser.add_argument(
        "--cutoff",
        type=float,
        default=-3,
        help="CLEAN cutoff (+ve = absolute, -ve = sigma).",
    )
    parser.add_argument(
        "--max_iter",
        type=int,
        default=10000,
        help="maximum number of CLEAN iterations.",
    )
    parser.add_argument("--gain", type=float, default=0.1, help="CLEAN loop gain.")
    parser.add_argument(
        "--window",
        type=float,
        default=None,
        help="Further CLEAN in mask to this threshold.",
    )

    return clean_parser


def cli():
    """Command-line interface"""

    from astropy.utils.exceptions import AstropyWarning

    warnings.simplefilter("ignore", category=AstropyWarning)
    from astropy.io.fits.verify import VerifyWarning

    warnings.simplefilter("ignore", category=VerifyWarning)

    gen_parser = generic_parser(parent_parser=True)
    common_parser = rmsynth_oncuts.rm_common_parser(parent_parser=True)
    rmclean_parser = clean_parser(parent_parser=True)

    parser = argparse.ArgumentParser(
        parents=[gen_parser, common_parser, rmclean_parser],
        formatter_class=UltimateHelpFormatter,
        description=rmclean_parser.description,
    )

    args = parser.parse_args()

    verbose = args.verbose
    rmv = args.rm_verbose
    host = args.host
    test_db(host=args.host, username=args.username, password=args.password)

    if rmv:
        logger.setLevel(
            level=logging.DEBUG,
        )
    elif verbose:
        logger.setLevel(
            level=logging.INFO,
        )
    main(
        field=args.field,
        outdir=Path(args.outdir),
        host=host,
        epoch=args.epoch,
        username=args.username,
        password=args.password,
        dimension=args.dimension,
        database=args.database,
        savePlots=args.save_plots,
        limit=args.limit,
        cutoff=args.cutoff,
        maxIter=args.max_iter,
        gain=args.gain,
        window=args.window,
        showPlots=args.show_plots,
        rm_verbose=args.rm_verbose,
    )


if __name__ == "__main__":
    cli()
