#!/usr/bin/env python3
"""Correct for the ionosphere in parallel"""

from __future__ import annotations

import argparse
import asyncio
import ftplib
import logging
import os
from pathlib import Path
from pprint import pformat
from typing import Any
from typing import NamedTuple as Struct

import astropy.units as u
import numpy as np
import pymongo
import requests
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time, TimeDelta
from prefect import flow, get_run_logger, task
from spinifex.download import download_or_copy_url
from spinifex.exceptions import IonexError
from spinifex.get_rm import get_rm_from_skycoord
from spinifex.image_tools import (
    correct_fits_images_from_table,
    get_integrated_rm,
    write_integrated_rm_table,
)
from spinifex.ionospheric.ionex_download import (
    new_cddis_format,
    old_cddis_format,
    should_use_new_cddis_format,
)
from spinifex.ionospheric.tec_data import CENTER_NAMES
from spinifex.times import get_unique_days
from tqdm.auto import tqdm

from arrakis.logger import TqdmToLogger, UltimateHelpFormatter, logger
from arrakis.utils.database import (
    get_db,
    get_field_db,
    get_field_observation,
    test_db,
    validate_sbid_field_pair,
)
from arrakis.utils.fitsutils import getfreq
from arrakis.utils.pipeline import generic_parser, logo_str, workdir_arg_parser

logger.setLevel(logging.INFO)
TQDM_OUT = TqdmToLogger(logger, level=logging.INFO)


class Prediction(Struct):
    """Spinifex prediction"""

    predict_file: str
    update: pymongo.UpdateOne


class FrionResults(Struct):
    prediction: Prediction
    correction: pymongo.UpdateOne


def parse_option_value(value: str) -> Any:
    """Turn a string with a unit (e.g. '30 min') into a Quantity, else leave it."""
    try:
        quantity = u.Quantity(value)
    except (TypeError, ValueError):
        return value
    return value if quantity.unit == u.dimensionless_unscaled else quantity


def parse_spinifex_options(options: list[str] | None) -> dict[str, Any]:
    """Parse KEY=VALUE strings into keyword arguments for the ionospheric model.

    Args:
        options (list[str] | None): KEY=VALUE strings

    Raises:
        ValueError: If an option is not KEY=VALUE

    Returns:
        dict[str, Any]: Keyword arguments
    """
    kwargs: dict[str, Any] = {}
    for option in options or []:
        key, sep, value = option.partition("=")
        if not sep:
            raise ValueError(f"Spinifex option '{option}' is not KEY=VALUE")
        kwargs[key.strip()] = parse_option_value(value.strip())
    return kwargs


def get_spinifex_server(server: str, formatter: str | None) -> tuple[str, str | None]:
    """Map FRion's IONEX server and formatter onto a Spinifex server and URL stem.

    The formatter picks the server layout, falling back to the server URL.
    'aiub' is not a Spinifex server, see `download_aiub_ionex`.

    Args:
        server (str): IONEX server URL
        formatter (str | None): FRion formatter name

    Raises:
        ValueError: If the layout is not recognised

    Returns:
        tuple[str, str | None]: Server name, and the URL stem if `server` matches it
    """
    layout = (formatter or server).lower()
    for keyword, name in (
        ("aiub", "aiub"),
        ("cddis", "cddis"),
        ("chapman", "chapman"),
        ("olsztyn", "igsiono"),
        ("igs", "igsiono"),
    ):
        if keyword in layout:
            url_stem = server.rstrip("/") if keyword in server.lower() else None
            return name, url_stem
    raise ValueError(f"Unknown IONEX server layout '{formatter or server}'")


def get_iono_kwargs(
    server: str,
    prefix: str,
    formatter: str | None,
    ionex_dir: Path,
    iono_model: str = "ionex",
    spinifex_options: list[str] | None = None,
) -> dict[str, Any]:
    """Build the Spinifex ionospheric model options.

    Args:
        server (str): IONEX server URL
        prefix (str): Spinifex analysis centre prefix (e.g. 'cod')
        formatter (str | None): FRion formatter name
        ionex_dir (Path): Directory to keep ionospheric data in
        iono_model (str, optional): Spinifex ionospheric model. Defaults to "ionex".
        spinifex_options (list[str] | None, optional): Extra KEY=VALUE options. Defaults to None.

    Returns:
        dict[str, Any]: Keyword arguments for `get_rm_from_skycoord`
    """
    iono_kwargs: dict[str, Any] = {"output_directory": ionex_dir}
    # TOMION options reject the IONEX server settings
    if iono_model.startswith("ionex"):
        name, url_stem = get_spinifex_server(server, formatter)
        # AIUB files are fetched by Arrakis and saved under CDDIS names
        iono_kwargs["server"] = "cddis" if name == "aiub" else name
        iono_kwargs["url_stem"] = url_stem if name != "aiub" else None
        iono_kwargs["prefix"] = prefix
    iono_kwargs.update(parse_spinifex_options(spinifex_options))
    return iono_kwargs


def download_aiub_ionex(times: Time, server: str, prefix: str, ionex_dir: Path) -> None:
    """Download IONEX files from the AIUB CODE archive.

    Files are saved under the names Spinifex expects from CDDIS, so Spinifex
    finds them on disk. The day after each observation day is included because
    Spinifex uses it to remove midnight jumps.

    Args:
        times (Time): Times of the prediction
        server (str): AIUB CODE directory URL
        prefix (str): Spinifex analysis centre prefix
        ionex_dir (Path): Directory to keep ionospheric data in
    """
    ionex_dir.mkdir(parents=True, exist_ok=True)
    for day in get_unique_days(Time([*times, *(times + 1 * u.day)])):
        if should_use_new_cddis_format(day, prefix):
            name = Path(new_cddis_format(day, prefix=prefix)).name
            aiub_name = name
        else:
            name = Path(old_cddis_format(day, prefix=prefix)).name
            aiub_name = name.upper()
        if (ionex_dir / name).exists():
            continue
        url = f"{server.rstrip('/')}/{day.datetime.year:04d}/{aiub_name}"
        downloaded = asyncio.run(download_or_copy_url(url, output_directory=ionex_dir))
        downloaded.rename(ionex_dir / name)


@task(name="FRion correction")
def correct_worker(
    beam: dict,
    outdir: Path,
    field: str,
    prediction: Prediction,
    island: dict,
) -> pymongo.UpdateOne:
    """Apply Spinifex corrections to a single island

    Args:
        beam (Dict): MongoDB beam document
        outdir (str): Output directory
        field (str): RACS field name
        prediction (Prediction): Spinifex prediction
        island (dict): RACS island document

    Returns:
        pymongo.UpdateOne: Pymongo update query
    """
    island_id = island["Source_ID"]
    qfile = outdir / str(beam["beams"][field]["q_file"])
    ufile = outdir / str(beam["beams"][field]["u_file"])

    qout = str(beam["beams"][field]["q_file"]).replace(".fits", ".ion.fits")
    uout = str(beam["beams"][field]["u_file"]).replace(".fits", ".ion.fits")

    # Writes <name>.ion.fits next to each input cube
    correct_fits_images_from_table(
        stokes_q_path=qfile,
        stokes_u_path=ufile,
        correction_path=Path(prediction.predict_file),
        suffix="ion",
        overwrite=True,
    )

    myquery = {"Source_ID": island_id}

    newvalues = {
        "$set": {
            f"beams.{field}.q_file_ion": qout,
            f"beams.{field}.u_file_ion": uout,
        }
    }
    return pymongo.UpdateOne(myquery, newvalues)


@task(name="FRion predction")
def predict_worker(
    island: dict,
    field: str,
    beam: dict,
    start_time: Time,
    end_time: Time,
    freq: np.ndarray,
    cutdir: Path,
    server: str = "ftp://ftp.aiub.unibe.ch/CODE/",
    prefix: str = "",
    formatter: str | None = None,
    proxy_server: str | None = None,
    pre_download: bool = False,
    iono_model: str = "ionex",
    spinifex_options: list[str] | None = None,
) -> Prediction:
    """Make a Spinifex prediction for a single island

    Args:
        island (Dict): Pymongo island document
        field (str): RACS field name
        beam (Dict): Pymongo beam document
        start_time (Time): Start time of the observation
        end_time (Time): End time of the observation
        freq (np.ndarray): Array of frequencies in Hz
        cutdir (str): Cutout directory
        server (str, optional): IONEX server URL
        prefix (str, optional): IONEX prefix, FRion style (e.g. 'codg')
        formatter (str | None, optional): IONEX server layout
        proxy_server (str | None, optional): Ignored, set HTTP(S)_PROXY instead
        pre_download (bool, optional): Ignored, files are always downloaded first
        iono_model (str, optional): Spinifex ionospheric model
        spinifex_options (list[str] | None, optional): Extra KEY=VALUE model options

    Returns:
        Prediction: Spinifex prediction file and pymongo update query
    """
    logger.setLevel(logging.INFO)

    if proxy_server is not None:
        logger.warning(
            "ionex_proxy_server is ignored, Spinifex uses the HTTP(S)_PROXY variables"
        )
    if pre_download:
        logger.info("IONEX files are always downloaded before the prediction")

    ifile: Path = cutdir / beam["beams"][field]["i_file"]
    i_dir = ifile.parent
    iname = island["Source_ID"]
    source = SkyCoord(ra=island["RA"] * u.deg, dec=island["Dec"] * u.deg)
    ionex_dir = cutdir.parent / "IONEXdata"
    # Hardcoded so compute nodes need no network. From the astropy sites registry.
    askap = EarthLocation.from_geodetic(
        lon=116.631424 * u.deg, lat=-26.697 * u.deg, height=360.99 * u.m
    )

    # At least 3 samples for the Simpson integration in Spinifex
    duration = (end_time - start_time).to(u.s)
    n_times = max(int(np.ceil(duration / (300 * u.s))) + 1, 3)
    times = start_time + np.linspace(0, 1, n_times) * duration

    # FRion style prefixes (e.g. 'codg') map onto Spinifex centres ('cod')
    prefixes_to_try = [
        _prefix[:3].lower()
        for _prefix in (prefix, "codg", "jplg", "casg", "esag", "upcg", "igsg")
    ]
    prefixes_to_try = [p for p in dict.fromkeys(prefixes_to_try) if p in CENTER_NAMES]
    for _prefix in prefixes_to_try:
        iono_kwargs = get_iono_kwargs(
            server=server,
            prefix=_prefix,
            formatter=formatter,
            ionex_dir=ionex_dir,
            iono_model=iono_model,
            spinifex_options=spinifex_options,
        )
        try:
            if iono_model.startswith("ionex") and (
                get_spinifex_server(server, formatter)[0] == "aiub"
            ):
                download_aiub_ionex(times, server, _prefix, ionex_dir)
            time_dep_rm = get_rm_from_skycoord(
                loc=askap,
                times=times,
                source=source,
                iono_model_name=iono_model,
                **iono_kwargs,
            )
            break
        except (*ftplib.all_errors, requests.RequestException, IonexError) as e:
            logger.error(f"Could not get IONEX data with prefix '{_prefix}': {e}")
            logger.warning("Trying next prefix.")
            continue
    else:
        raise FileNotFoundError(
            f"Could not find IONEX file with prefixes {prefixes_to_try}"
        )

    freq_hz = np.asarray(freq) * u.Hz
    integrated_rm = get_integrated_rm(time_dep_rm=time_dep_rm, freq_arr=freq_hz)

    predict_file = os.path.join(i_dir, f"{iname}_ion.fits")
    write_integrated_rm_table(
        integrated_rm, freq_hz, Path(predict_file), overwrite=True
    )
    logger.info(f"Prediction file: {predict_file}")

    myquery = {"Source_ID": iname}

    newvalues = {
        "$set": {
            "frion": {
                "times": (
                    time_dep_rm.times.mjd * 86400.0
                ).tolist(),  # Turn back into MJD seconds for backwards compatibility
                "RMs": np.asarray(time_dep_rm.rm).tolist(),
                "theta_real": integrated_rm.theta.real.tolist(),
                "theta_imag": integrated_rm.theta.imag.tolist(),
            }
        }
    }

    update = pymongo.UpdateOne(myquery, newvalues)

    return Prediction(predict_file, update)


@task(name="Index beams")
def index_beams(island: dict, beams: list[dict]) -> dict:
    island_id = island["Source_ID"]
    beam_idx = [i for i, b in enumerate(beams) if b["Source_ID"] == island_id][0]
    beam = beams[beam_idx]
    return beam


# We reduce the inner loop to a serial call
# This is to avoid overwhelming the Prefect server
@task(name="FRion loop")
def serial_loop(
    island: dict,
    field: str,
    beam: dict,
    start_time: Time,
    end_time: Time,
    freq_hz_array: np.ndarray,
    cutdir: Path,
    ionex_server: str,
    ionex_prefix: str,
    ionex_proxy_server: str | None,
    ionex_formatter: str | None,
    ionex_predownload: bool,
    iono_model: str = "ionex",
    spinifex_options: list[str] | None = None,
) -> FrionResults:
    prediction = predict_worker.fn(
        island=island,
        field=field,
        beam=beam,
        start_time=start_time,
        end_time=end_time,
        freq=freq_hz_array,
        cutdir=cutdir,
        server=ionex_server,
        prefix=ionex_prefix,
        proxy_server=ionex_proxy_server,
        formatter=ionex_formatter,
        pre_download=ionex_predownload,
        iono_model=iono_model,
        spinifex_options=spinifex_options,
    )
    correction = correct_worker.fn(
        beam=beam,
        outdir=cutdir,
        field=field,
        prediction=prediction,
        island=island,
    )

    return FrionResults(prediction=prediction, correction=correction)


@flow(name="FRion")
def main(
    field: str,
    outdir: Path,
    host: str,
    epoch: int,
    sbid: int | None = None,
    username: str | None = None,
    password: str | None = None,
    database: bool = False,
    ionex_server: str = "ftp://ftp.aiub.unibe.ch/CODE/",
    ionex_prefix: str = "codg",
    ionex_proxy_server: str | None = None,
    ionex_formatter: str | None = "ftp.aiub.unibe.ch",
    ionex_predownload: bool = False,
    limit: int | None = None,
    iono_model: str = "ionex",
    spinifex_options: list[str] | None = None,
) -> None:
    """Ionospheric Faraday rotation flow, using Spinifex

    Args:
        field (str): RACS field name
        outdir (Path): Output directory
        host (str): MongoDB host IP address
        epoch (int): Epoch of observation
        sbid (int, optional): SBID of observation. Defaults to None.
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo password. Defaults to None.
        database (bool, optional): Update database. Defaults to False.
        ionex_server (str, optional): IONEX server. Defaults to "ftp://ftp.aiub.unibe.ch/CODE/".
        ionex_prefix (str, optional): IONEX prefix. Defaults to "codg".
        ionex_proxy_server (str, optional): Ignored, set HTTP(S)_PROXY instead. Defaults to None.
        ionex_formatter (str, optional): IONEX server layout. Defaults to "ftp.aiub.unibe.ch".
        ionex_predownload (bool, optional): Ignored, files are always downloaded first. Defaults to False.
        limit (int, optional): Limit to number of islands. Defaults to None.
        iono_model (str, optional): Spinifex ionospheric model. Defaults to "ionex".
        spinifex_options (list[str], optional): Extra KEY=VALUE model options. Defaults to None.
    """
    logger = get_run_logger()
    # Query database for data
    outdir = outdir.absolute()
    cutdir = outdir / "cutouts"

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
            raise ValueError(f"SBID {sbid} does not match field {field}")

    query_1 = {"$and": [{f"beams.{field}": {"$exists": True}}]}

    if sbid is not None:
        query_1["$and"].append({f"beams.{field}.SBIDs": sbid})

    beams = list(beams_col.find(query_1).sort("Source_ID"))
    island_ids = sorted(beams_col.distinct("Source_ID", query_1))

    # Get FRion arguments
    query_2 = {"Source_ID": {"$in": island_ids}}
    islands = list(island_col.find(query_2).sort("Source_ID"))

    field_col = get_field_db(
        host=host, epoch=epoch, username=username, password=password
    )
    # Use the processed SBID's observation if given, otherwise the field's
    # SELECT=1 observation, so a re-observed field gets the right times
    field_data = get_field_observation(field_col, field_name=field, sbid=sbid)
    logger.info(f"{field_data=}")

    start_time = Time(field_data["SCAN_START"] * u.second, format="mjd")
    end_time = start_time + TimeDelta(field_data["SCAN_TINT"] * u.second)

    freq = getfreq(
        os.path.join(cutdir, f"{beams[0]['beams'][f'{field}']['q_file']}"),
    )

    if limit is not None:
        logger.info(f"Limiting to {limit} islands")
        islands = islands[:limit]

    beams_cor = []
    for island in islands:
        island_id = island["Source_ID"]
        beam_idx = [i for i, b in enumerate(beams) if b["Source_ID"] == island_id][0]
        beam = beams[beam_idx]
        beams_cor.append(beam)

    # do one prediction to get the IONEX files
    _ = predict_worker(
        island=islands[0],
        field=field,
        beam=beams_cor[0],
        start_time=start_time,
        end_time=end_time,
        freq=freq.to(u.Hz).value,
        cutdir=cutdir,
        server=ionex_server,
        prefix=ionex_prefix,
        proxy_server=ionex_proxy_server,
        formatter=ionex_formatter,
        pre_download=ionex_predownload,
        iono_model=iono_model,
        spinifex_options=spinifex_options,
    )

    frion_results = []
    assert len(islands) == len(beams_cor), "Islands and beams must be the same length"
    for island, beam in tqdm(
        zip(islands, beams_cor),
        desc="Submitting tasks",
        file=TQDM_OUT,
        total=len(islands),
    ):
        frion_result = serial_loop.submit(
            island=island,
            field=field,
            beam=beam,
            start_time=start_time,
            end_time=end_time,
            freq_hz_array=freq.to(u.Hz).value,
            cutdir=cutdir,
            ionex_server=ionex_server,
            ionex_prefix=ionex_prefix,
            ionex_proxy_server=ionex_proxy_server,
            ionex_formatter=ionex_formatter,
            ionex_predownload=ionex_predownload,
            iono_model=iono_model,
            spinifex_options=spinifex_options,
        )
        frion_results.append(frion_result)

    predictions = []
    corrections = []
    for result in frion_results:
        predictions.append(result.result().prediction)
        corrections.append(result.result().correction)

    updates_arrays = [p.update for p in predictions]
    updates = corrections
    if database:
        logger.info("Updating beams database...")
        db_res = beams_col.bulk_write(updates, ordered=False)
        logger.info(pformat(db_res.bulk_api_result))

        logger.info("Updating island database...")
        db_res = island_col.bulk_write(updates_arrays, ordered=False)
        logger.info(pformat(db_res.bulk_api_result))


def frion_parser(parent_parser: bool = False) -> argparse.ArgumentParser:
    # Help string to be shown using the -h option
    descStr = f"""
    {logo_str}
    Arrakis Stage:
    Correct for ionospheric Faraday rotation using Spinifex

    """

    # Parse the command line options
    frion_parser = argparse.ArgumentParser(
        add_help=not parent_parser,
        description=descStr,
        formatter_class=UltimateHelpFormatter,
    )
    parser = frion_parser.add_argument_group("frion arguments")

    parser.add_argument(
        "--ionex_server",
        type=str,
        default="ftp://ftp.aiub.unibe.ch/CODE/",
        help="IONEX server",
    )

    parser.add_argument(
        "--ionex_prefix",
        type=str,
        default="codg",
    )

    parser.add_argument(
        "--ionex_formatter",
        type=str,
        default="ftp.aiub.unibe.ch",
        help="IONEX formatter.",
    )

    parser.add_argument(
        "--ionex_proxy_server",
        type=str,
        default=None,
        help="Proxy server. Ignored, Spinifex uses the HTTP(S)_PROXY variables.",
    )

    parser.add_argument(
        "--ionex_predownload",
        action="store_true",
        help="Pre-download IONEX files.",
    )

    parser.add_argument(
        "--iono_model",
        type=str,
        default="ionex",
        choices=["ionex", "ionex_iri", "tomion", "tomion_dual"],
        help="Spinifex ionospheric model.",
    )

    parser.add_argument(
        "--spinifex_options",
        type=str,
        nargs="+",
        default=None,
        metavar="KEY=VALUE",
        help="Extra Spinifex ionospheric model options, e.g. 'solution=rapid' "
        "'height=350 km'. These override the ionex arguments.",
    )

    return frion_parser


def cli():
    """Command-line interface"""
    import warnings

    from astropy.utils.exceptions import AstropyWarning

    warnings.simplefilter("ignore", category=AstropyWarning)
    from astropy.io.fits.verify import VerifyWarning

    warnings.simplefilter("ignore", category=VerifyWarning)
    warnings.simplefilter("ignore", category=RuntimeWarning)

    gen_parser = generic_parser(parent_parser=True)
    work_parser = workdir_arg_parser(parent_parser=True)
    f_parser = frion_parser(parent_parser=True)
    parser = argparse.ArgumentParser(
        parents=[gen_parser, work_parser, f_parser],
        formatter_class=UltimateHelpFormatter,
        description=f_parser.description,
    )
    args = parser.parse_args()

    verbose = args.verbose
    if verbose:
        logger.setLevel(logging.INFO)

    test_db(host=args.host, username=args.username, password=args.password)

    main(
        field=args.field,
        sbid=args.sbid,
        outdir=Path(
            args.datadir,
        ),
        host=args.host,
        epoch=args.epoch,
        username=args.username,
        password=args.password,
        database=args.database,
        ionex_server=args.ionex_server,
        ionex_proxy_server=args.ionex_proxy_server,
        ionex_formatter=args.ionex_formatter,
        ionex_prefix=args.ionex_prefix,
        ionex_predownload=args.ionex_predownload,
        iono_model=args.iono_model,
        spinifex_options=args.spinifex_options,
    )


if __name__ == "__main__":
    cli()
