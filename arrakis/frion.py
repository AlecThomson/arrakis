#!/usr/bin/env python3
"""Correct for the ionosphere in parallel"""
import logging
import os
import time
from glob import glob
from pprint import pformat
from shutil import copyfile
from typing import Dict, List, Optional, Tuple, Union

import astropy.units as u
import dask
import numpy as np
import pymongo
from astropy.time import Time, TimeDelta
from dask import delayed
from dask.distributed import Client, LocalCluster, progress, wait
from FRion import correct, predict

from arrakis.logger import logger
from arrakis.utils import get_db, get_field_db, getfreq, test_db, tqdm_dask, try_mkdir

logger.setLevel(logging.INFO)


@delayed
def correct_worker(
    beam: Dict, outdir: str, field: str, predict_file: str, island_id: str
) -> pymongo.UpdateOne:
    """Apply FRion corrections to a single island

    Args:
        beam (Dict): MongoDB beam document
        outdir (str): Output directory
        field (str): RACS field name
        predict_file (str): FRion prediction file
        island_id (str): RACS island ID

    Returns:
        pymongo.UpdateOne: Pymongo update query
    """
    qfile = os.path.join(outdir, beam["beams"][field]["q_file"])
    ufile = os.path.join(outdir, beam["beams"][field]["u_file"])

    qout = beam["beams"][field]["q_file"].replace(".fits", ".ion.fits")
    uout = beam["beams"][field]["u_file"].replace(".fits", ".ion.fits")

    qout_f = os.path.join(outdir, qout)
    uout_f = os.path.join(outdir, uout)

    correct.apply_correction_to_files(
        qfile, ufile, predict_file, qout_f, uout_f, overwrite=True
    )

    myquery = {"Source_ID": island_id}

    newvalues = {
        "$set": {
            f"beams.{field}.q_file_ion": qout,
            f"beams.{field}.u_file_ion": uout,
        }
    }
    return pymongo.UpdateOne(myquery, newvalues)


@delayed(nout=2)
def predict_worker(
    island: Dict,
    field: str,
    beam: Dict,
    start_time: Time,
    end_time: Time,
    freq: np.ndarray,
    cutdir: str,
    plotdir: str,
    server: str = "ftp://ftp.aiub.unibe.ch/CODE/",
    proxy_server: Optional[str] = None,
) -> Tuple[str, pymongo.UpdateOne]:
    """Make FRion prediction for a single island

    Args:
        island (Dict): Pymongo island document
        field (str): RACS field name
        beam (Dict): Pymongo beam document
        start_time (Time): Start time of the observation
        end_time (Time): End time of the observation
        freq (np.ndarray): Array of frequencies with units
        cutdir (str): Cutout directory
        plotdir (str): Plot directory

    Returns:
        Tuple[str, pymongo.UpdateOne]: FRion prediction file and pymongo update query
    """
    logger.setLevel(logging.INFO)

    ifile = os.path.join(cutdir, beam["beams"][field]["i_file"])
    i_dir = os.path.dirname(ifile)
    iname = island["Source_ID"]
    ra = island["RA"]
    dec = island["Dec"]

    # Tricking the ionex lookup to use the a custom server
    proxy_args = {
        "proxy_type": None,
        "proxy_port": None,
        "proxy_user": None,
        "proxy_pass": None,
    }
    logger.info(f"Set up empty Proxy structure.")

    times, RMs, theta = predict.calculate_modulation(
        start_time=start_time.fits,
        end_time=end_time.fits,
        freq_array=freq,
        telescope_location=predict.get_telescope_coordinates("ASKAP"),
        ra=ra,
        dec=dec,
        timestep=300.0,
        ionexPath=os.path.join(os.path.dirname(cutdir), "IONEXdata"),
        server=server,
        proxy_server=proxy_server,
        **proxy_args,
    )
    predict_file = os.path.join(i_dir, f"{iname}_ion.txt")
    predict.write_modulation(freq_array=freq, theta=theta, filename=predict_file)

    plot_file = os.path.join(i_dir, f"{iname}_ion.pdf")
    predict.generate_plots(
        times, RMs, theta, freq, position=[ra, dec], savename=plot_file
    )
    plot_files = glob(os.path.join(i_dir, "*ion.pdf"))
    logger.info(f"Plotting files: {plot_files=}")
    for src in plot_files:
        base = os.path.basename(src)
        dst = os.path.join(plotdir, base)
        copyfile(src, dst)

    myquery = {"Source_ID": iname}

    newvalues = {
        "$set": {
            "frion": {
                "times": (
                    times.mjd * 86400.0
                ).tolist(),  # Turn back into MJD seconds for backwards compatibility
                "RMs": RMs.tolist(),
                "theta_real": theta.real.tolist(),
                "theta_imag": theta.imag.tolist(),
            }
        }
    }

    update = pymongo.UpdateOne(myquery, newvalues)

    return predict_file, update


def main(
    field: str,
    outdir: str,
    host: str,
    username: Optional[str] = None,
    password: Optional[str] = None,
    database=False,
    verbose=True,
    ionex_server: str = "ftp://ftp.aiub.unibe.ch/CODE/",
    ionex_proxy_server: Optional[str] = None,
):
    """Main script

    Args:
        field (str): RACS field name
        outdir (str): Output directory
        host (str): MongoDB host IP address
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo passwrod. Defaults to None.
        database (bool, optional): Update database. Defaults to False.
        verbose (bool, optional): Verbose output. Defaults to True.
        ionex_server (str, optional): IONEX server. Defaults to "ftp://ftp.aiub.unibe.ch/CODE/".
        ionex_proxy_server (str, optional): Proxy server. Defaults to None.
    """
    # Query database for data
    outdir = os.path.abspath(outdir)
    cutdir = os.path.join(outdir, "cutouts")

    plotdir = os.path.join(cutdir, "plots")
    try_mkdir(plotdir)

    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )

    query_1 = {"$and": [{f"beams.{field}": {"$exists": True}}]}

    beams = list(beams_col.find(query_1).sort("Source_ID"))
    island_ids = sorted(beams_col.distinct("Source_ID", query_1))

    # Get FRion arguments
    query_2 = {"Source_ID": {"$in": island_ids}}
    islands = list(island_col.find(query_2).sort("Source_ID"))

    field_col = get_field_db(host, username=username, password=password)
    query_3 = {"FIELD_NAME": f"{field}"}
    logger.info(f"{query_3}")

    # Get most recent SBID
    if field_col.count_documents(query_3) > 1:
        field_datas = list(field_col.find({"FIELD_NAME": f"{field}"}))
        sbids = [f["CAL_SBID"] for f in field_datas]
        max_idx = np.argmax(sbids)
        logger.info(f"Using CAL_SBID {sbids[max_idx]}")
        field_data = field_datas[max_idx]
    else:
        field_data = field_col.find_one({"FIELD_NAME": f"{field}"})

    logger.info(f"{field_data=}")

    start_time = Time(field_data["SCAN_START"] * u.second, format="mjd")
    end_time = start_time + TimeDelta(field_data["SCAN_TINT"] * u.second)

    freq = getfreq(
        os.path.join(cutdir, f"{beams[0]['beams'][f'{field}']['q_file']}"),
    )  # Type: u.Quantity

    # Loop over islands in parallel
    outputs = []
    updates_arrays = []
    for island in islands:
        island_id = island["Source_ID"]
        beam_idx = [i for i, b in enumerate(beams) if b["Source_ID"] == island_id][0]
        beam = beams[beam_idx]
        # Get FRion predictions
        predict_file, update = predict_worker(
            island=island,
            field=field,
            beam=beam,
            start_time=start_time,
            end_time=end_time,
            freq=freq.to(u.Hz).value,
            cutdir=cutdir,
            plotdir=plotdir,
            server=ionex_server,
            proxy_server=ionex_proxy_server,
        )
        updates_arrays.append(update)
        # Apply FRion predictions
        output = correct_worker(
            beam=beam,
            outdir=cutdir,
            field=field,
            predict_file=predict_file,
            island_id=island_id,
        )
        outputs.append(output)
    # Wait for IONEX data I guess...
    _ = outputs[0].compute()
    time.sleep(10)
    # Execute
    futures, future_arrays = dask.persist(outputs, updates_arrays)
    # dumb solution for https://github.com/dask/distributed/issues/4831
    time.sleep(10)
    tqdm_dask(
        futures, desc="Running FRion", disable=(not verbose), total=len(islands) * 3
    )
    if database:
        logger.info("Updating beams database...")
        updates = [f.compute() for f in futures]
        db_res = beams_col.bulk_write(updates, ordered=False)
        logger.info(pformat(db_res.bulk_api_result))

        logger.info("Updating island database...")
        updates_arrays_cmp = [f.compute() for f in future_arrays]

        db_res = island_col.bulk_write(updates_arrays_cmp, ordered=False)
        logger.info(pformat(db_res.bulk_api_result))


def cli():
    """Command-line interface"""
    import argparse
    import warnings

    from astropy.utils.exceptions import AstropyWarning

    warnings.simplefilter("ignore", category=AstropyWarning)
    from astropy.io.fits.verify import VerifyWarning

    warnings.simplefilter("ignore", category=VerifyWarning)
    warnings.simplefilter("ignore", category=RuntimeWarning)
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
    Arrakis Stage:
    Correct for ionospheric Faraday rotation

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field to mosaic - e.g. 2132-50A."
    )
    parser.add_argument(
        "outdir",
        metavar="outdir",
        type=str,
        help="Directory containing cutouts (in subdir outdir/cutouts).",
    )

    parser.add_argument(
        "host",
        metavar="host",
        type=str,
        help="Host of mongodb (probably $hostname -i).",
    )

    parser.add_argument(
        "--username", type=str, default="admin", help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )

    parser.add_argument(
        "-m", "--database", action="store_true", help="Add data to MongoDB [False]."
    )

    parser.add_argument(
        "-s",
        "--ionex_server",
        type=str,
        default="ftp://ftp.aiub.unibe.ch/CODE/",
        help="IONEX server [ftp://ftp.aiub.unibe.ch/CODE/].",
    )

    parser.add_argument(
        "-p",
        "--ionex_proxy_server",
        type=str,
        default=None,
        help="Proxy server [None].",
    )

    parser.add_argument(
        "-v", dest="verbose", action="store_true", help="verbose output [False]."
    )

    args = parser.parse_args()

    verbose = args.verbose
    if verbose:
        logger.setLevel(logging.INFO)

    cluster = LocalCluster(
        n_workers=10, processes=True, threads_per_worker=1, local_directory="/dev/shm"
    )
    client = Client(cluster)
    logger.info(client)

    test_db(
        host=args.host, username=args.username, password=args.password, verbose=verbose
    )

    main(
        field=args.field,
        outdir=args.outdir,
        host=args.host,
        username=args.username,
        password=args.password,
        database=args.database,
        verbose=verbose,
        ionex_server=args.ionex_server,
        ionex_proxy_server=args.ionex_proxy_server,
    )


if __name__ == "__main__":
    cli()
