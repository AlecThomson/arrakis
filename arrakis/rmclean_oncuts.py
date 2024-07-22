#!/usr/bin/env python3
"""Run RM-synthesis on cutouts in parallel"""

from __future__ import annotations

import argparse
import logging
import os
import warnings
from pathlib import Path
from pprint import pformat
from shutil import copyfile

import matplotlib as mpl
import numpy as np
import pymongo
from matplotlib import pyplot as plt
from prefect import flow, task
from RMtools_1D import do_RMclean_1D
from RMtools_3D import do_RMclean_3D
from tqdm.auto import tqdm

from arrakis import rmsynth_oncuts
from arrakis.logger import TqdmToLogger, UltimateHelpFormatter, logger
from arrakis.utils.database import (
    get_db,
    get_field_db,
    test_db,
    validate_sbid_field_pair,
)
from arrakis.utils.pipeline import generic_parser, logo_str, workdir_arg_parser

mpl.use("Agg")
logger.setLevel(logging.INFO)
TQDM_OUT = TqdmToLogger(logger, level=logging.INFO)


@task(name="1D RM-CLEAN")
def rmclean1d(
    field: str,
    comp: dict,
    outdir: Path,
    cutoff: float = -3,
    maxIter=10000,
    gain=0.1,
    sbid: int | None = None,
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
        savePlots (bool, optional): Save CLEAN plots. Defaults to False.
        rm_verbose (bool, optional): Verbose RM-CLEAN. Defaults to True.

    Returns:
        pymongo.UpdateOne: MongoDB update query.
    """
    cname = comp["Gaussian_ID"]
    logger.debug(f"Working on {comp}")
    save_name = field if sbid is None else f"{field}_{sbid}"

    assert (
        comp["rm_outputs_1d"]["field"] == save_name
    ), f"Field mismatch - expected {save_name} got {comp['rm_outputs_1d']['field']}"

    rm1dfiles = comp["rm_outputs_1d"]["rm1dfiles"]
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
            msg = f"File does not exist: '{f}'"
            raise FileNotFoundError(msg)

    nBits = 32
    try:
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
            showPlots=False,
            verbose=rm_verbose,
            prefixOut=prefix,
            saveFigures=False,
            window=window,
        )
    except Exception as e:
        logger.error(f"Error running RM-CLEAN on {cname}: {e}")
        myquery = {"Gaussian_ID": cname, "rm_outputs_1d.field": save_name}
        operation = {"$set": {"rm_outputs_1d.$.rmclean1d": False}}
        return pymongo.UpdateOne(myquery, operation, upsert=True)

    # Ensure JSON serializable
    for k, v in outdict.items():
        if isinstance(v, (np.float64, np.float32)):
            outdict[k] = float(v)
        elif isinstance(v, (np.int_, np.int32)):
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
    myquery = {
        "Gaussian_ID": cname,
        "rm_outputs_1d.field": save_name,
    }

    to_update = {
        "rm_outputs_1d.$.rmclean1d": True,
        "rm_outputs_1d.$.rmclean_summary": outdict,
    }
    newvalues = {"$set": to_update}
    return pymongo.UpdateOne(myquery, newvalues, upsert=True)


@task(name="3D RM-CLEAN")
def rmclean3d(
    field: str,
    island: dict,
    outdir: Path,
    sbid: int | None = None,
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
    rm3dfiles = island["rm_outputs_3d"]["rm3dfiles"]
    save_name = field if sbid is None else f"{field}_{sbid}"
    assert (
        island["rm_outputs_3d"]["field"] == save_name
    ), f"Field mismatch - expected {save_name} got {island['rm_outputs_3d']['field']}"

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
    myquery = {"Source_ID": iname, "rm_outputs_3d.field": save_name}
    to_update = {
        "rm_outputs_3d.$.rmclean3d": True,
    }
    newvalues = {"$set": to_update}

    return pymongo.UpdateOne(myquery, newvalues, upsert=True)


@flow(name="RM-CLEAN on cutouts")
def main(
    field: str,
    outdir: Path,
    host: str,
    epoch: int,
    sbid: int | None = None,
    username: str | None = None,
    password: str | None = None,
    dimension="1d",
    database=False,
    savePlots=True,
    limit: int | None = None,
    cutoff: float = -3,
    maxIter=10000,
    gain=0.1,
    window=None,
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
        rm_verbose (bool, optional): Verbose output from RM-CLEAN. Defaults to False.
    """
    outdir = outdir.absolute() / "cutouts"

    # default connection (ie, local)
    beams_col, island_col, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )

    # Check for SBID match
    if sbid is not None:
        field_col = get_field_db(
            host=host,
            epoch=epoch,
            username=username,
            password=password,
        )
        sbid_check = validate_sbid_field_pair(
            field_name=field,
            sbid=sbid,
            field_col=field_col,
        )
        if not sbid_check:
            msg = f"SBID {sbid} does not match field {field}"
            raise ValueError(msg)

    query = {"$and": [{f"beams.{field}": {"$exists": True}}]}
    if sbid is not None:
        query["$and"].append({f"beams.{field}.SBIDs": sbid})

    all_island_ids = sorted(beams_col.distinct("Source_ID", query))

    save_name = field if sbid is None else f"{field}_{sbid}"
    if dimension == "3d":
        query = {
            "$and": [
                {"Source_ID": {"$in": all_island_ids}},
                {"rm_outputs_3d.field": save_name, "rm_outputs_3d.rmsynth3d": True},
            ]
        }
        # exit()
        pipeline = [
            {"$match": query},
            {
                "$project": {
                    "Source_ID": 1,
                    "Gaussian_ID": 1,
                    "rm_outputs_3d": {
                        "$arrayElemAt": [
                            "$rm_outputs_3d",
                            {"$indexOfArray": ["$rm_outputs_3d.field", save_name]},
                        ]
                    },
                }
            },
        ]
        islands = list(island_col.aggregate(pipeline))
        n_comp = island_col.count_documents(query)
        update_3d = {"rm_outputs_3d.$.rmclean3d": False}
        operation_3d = {"$set": update_3d}
        logger.info(pformat(operation_3d))

        result = island_col.update_many(
            query,
            operation_3d,
        )
        logger.info(pformat(result.raw_result))
        if limit is not None:
            count = limit
            n_island = count
            islands = islands[:count]

    elif dimension == "1d":
        query = {
            "$and": [
                {"Source_ID": {"$in": all_island_ids}},
                {
                    "rm_outputs_1d": {
                        "$elemMatch": {"field": save_name, "rmsynth1d": True}
                    }
                },
                # {"rm_outputs_1d.field": save_name, "rm_outputs_1d.rmsynth1d": True},
            ]
        }
        # exit()
        pipeline = [
            {"$match": query},
            {
                "$project": {
                    "Source_ID": 1,
                    "Gaussian_ID": 1,
                    "rm_outputs_1d": {
                        "$arrayElemAt": [
                            "$rm_outputs_1d",
                            {"$indexOfArray": ["$rm_outputs_1d.field", save_name]},
                        ]
                    },
                }
            },
        ]
        components = list(comp_col.aggregate(pipeline))
        n_comp = comp_col.count_documents(query)
        update_1d = {"rm_outputs_1d.$.rmclean1d": False}
        # filter_condition = [{"elem.field": save_name}]
        operation_1d = {"$set": update_1d}
        logger.info(pformat(operation_1d))

        result = comp_col.update_many(
            query,
            operation_1d,
            upsert=True,
        )
        logger.info(pformat(result.raw_result))

        if limit is not None:
            count = limit
            n_comp = count
            components = components[:count]

    if dimension == "1d":
        outputs = []
        logger.info(f"Running RM-CLEAN on {n_comp} components")
        for comp in tqdm(components, total=n_comp, desc="RM-CLEAN 1D", file=TQDM_OUT):
            output = rmclean1d.submit(
                comp=comp,
                field=field,
                sbid=sbid,
                outdir=outdir,
                cutoff=cutoff,
                maxIter=maxIter,
                gain=gain,
                savePlots=savePlots,
                rm_verbose=rm_verbose,
                window=window,
            )
            outputs.append(output)

    elif dimension == "3d":
        logger.info(f"Running RM-CLEAN on {n_island} islands")
        outputs = []
        for island in tqdm(islands, total=n_island, desc="RM-CLEAN 3D", file=TQDM_OUT):
            output = rmclean3d.submit(
                field=field,
                island=island,
                sbid=sbid,
                outdir=outdir,
                cutoff=cutoff,
                maxIter=maxIter,
                gain=gain,
                rm_verbose=rm_verbose,
            )
            outputs.append(output)

    else:
        msg = f"Dimension {dimension} not supported."
        raise ValueError(msg)

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
    work_parser = workdir_arg_parser(parent_parser=True)
    common_parser = rmsynth_oncuts.rm_common_parser(parent_parser=True)
    rmclean_parser = clean_parser(parent_parser=True)

    parser = argparse.ArgumentParser(
        parents=[gen_parser, work_parser, common_parser, rmclean_parser],
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
        sbid=args.sbid,
        outdir=Path(args.datadir),
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
        rm_verbose=args.rm_verbose,
    )


if __name__ == "__main__":
    cli()
