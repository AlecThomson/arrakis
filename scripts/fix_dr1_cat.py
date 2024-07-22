#!/usr/bin/env python3
"""Post process DR1 catalog"""

import logging
import os
import pickle
from pathlib import Path

import astropy.units as u
import numpy as np
from arrakis.logger import logger
from arrakis.makecat import (
    compute_local_rm_flag,
    get_fit_func,
    is_leakage,
    write_votable,
)
from astropy.coordinates import SkyCoord
from astropy.table import Column, Table
from astropy.time import Time
from astropy.units import cds
from rmtable import RMTable

from spica import SPICA


def fix_fields(
    tab: Table,
    survey_dir: Path,
    epoch: int = 0,
) -> Table:
    # Get field data, and index by field/tile ID
    field_path = survey_dir / "db" / f"epoch_{epoch}" / "field_data.csv"
    field = Table.read(field_path)
    field = field[field["SELECT"] == 1]
    field.add_index("FIELD_NAME")
    tab.add_index("tile_id")

    # Compare the fields we have to those we want
    fields_in_cat = list(set(tab["tile_id"]))
    fields_in_spica = [f"RACS_{name}" for name in SPICA]
    logger.debug(f"Fields in catalogue: {fields_in_cat}")
    logger.debug(f"Fields in spica: {fields_in_spica}")
    fields_not_in_spica = [f for f in fields_in_cat if f not in fields_in_spica]
    spica_field = field.loc[fields_in_spica]
    spica_field_coords = SkyCoord(
        spica_field["RA_DEG"], spica_field["DEC_DEG"], unit=(u.deg, u.deg), frame="icrs"
    )
    start_times = Time(spica_field["SCAN_START"] * u.second, format="mjd")
    spica_field.add_column(
        Column(
            start_times.to_value("mjd"),
            name="start_time",
            unit=cds.MJD,
        ),
    )
    # These are the sources to update
    sources_to_fix = tab.loc[fields_not_in_spica]
    logger.info(f"Found {len(sources_to_fix)} sources to fix")

    source_coords = SkyCoord(sources_to_fix["ra"], sources_to_fix["dec"])

    # Get separation between source and field centres
    seps = []
    for c in spica_field_coords:
        sep = c.separation(source_coords)
        seps.append(sep.to(u.deg).value)
    # Find the closest field and set the tile_id etc in catalogue
    sep_arr = np.array(seps) * u.deg
    min_idx = np.argmin(sep_arr, axis=0)
    min_seps = np.min(sep_arr, axis=0)
    closest_fields = np.array(fields_in_spica)[min_idx]
    new_tab = tab.copy()
    idx = new_tab.loc_indices[fields_not_in_spica]

    # Update tile_id, SBID, start time, and field sep
    new_tab.remove_indices("tile_id")

    all_fields = new_tab["tile_id"].value
    all_fields[idx] = closest_fields
    new_tab.replace_column(
        "tile_id",
        Column(
            all_fields,
            name="tile_id",
        ),
    )

    all_seps = (
        new_tab["separation_tile_centre"].value * new_tab["separation_tile_centre"].unit
    )
    all_seps[idx] = min_seps

    all_sbids = new_tab["sbid"].value
    all_sbids[idx] = spica_field["SBID"][min_idx].value

    all_start_times = new_tab["start_time"]
    all_start_times[idx] = spica_field["start_time"][min_idx]

    # Update the columns
    new_tab.replace_column(
        "separation_tile_centre",
        Column(
            data=all_seps,
            name="separation_tile_centre",
            unit=all_seps.unit,
        ),
    )
    new_tab.replace_column(
        "beamdist",
        Column(
            data=all_seps,
            name="beamdist",
            unit=all_seps.unit,
        ),
    )

    new_tab.replace_column(
        "sbid",
        Column(
            data=all_sbids,
            name="sbid",
        ),
    )

    new_tab.replace_column(
        "start_time",
        Column(
            data=all_start_times,
            name="start_time",
            unit=all_start_times.unit,
        ),
    )

    # Fix the units - Why does VOTable do this?? Thanks I hate it
    dumb_units = {
        "Jy.beam-1": u.Jy / u.beam,
        "mJy.beam-1": u.mJy / u.beam,
        "day": u.d,
    }
    for col in new_tab.colnames:
        if str(new_tab[col].unit) in dumb_units.keys():
            new_unit = dumb_units[str(new_tab[col].unit)]
            logger.debug(f"Fixing {col} unit from {new_tab[col].unit} to {new_unit}")
            new_tab[col].unit = new_unit
            new_tab.units[col] = new_unit

    # Convert all mJy to Jy
    for col in new_tab.colnames:
        if new_tab[col].unit == u.mJy:
            logger.debug(f"Converting {col} unit from {new_tab[col].unit} to {u.Jy}")
            new_tab[col] = new_tab[col].to(u.Jy)
            new_tab.units[col] = u.Jy
        if new_tab[col].unit == u.mJy / u.beam:
            logger.debug(
                f"Converting {col} unit from {new_tab[col].unit} to {u.Jy / u.beam}"
            )
            new_tab[col] = new_tab[col].to(u.Jy / u.beam)
            new_tab.units[col] = u.Jy / u.beam

    return new_tab


def main(cat: str, survey_dir: Path, epoch: int = 0):
    logger.info(f"Reading {cat}")
    tab = RMTable.read(cat)
    logger.info(f"Fixing {cat}")

    fix_tab = fix_fields(tab=tab, survey_dir=survey_dir, epoch=epoch)
    fit, fig = get_fit_func(fix_tab, do_plot=True, nbins=16, degree=4)
    fig.savefig("leakage_fit_dr1_fix.pdf")
    leakage_flag = is_leakage(
        fix_tab["fracpol"].value, fix_tab["beamdist"].to(u.deg).value, fit
    )
    fix_tab.replace_column(
        "leakage_flag",
        Column(
            leakage_flag,
            name="leakage_flag",
        ),
    )
    leakage = fit(fix_tab["separation_tile_centre"].to(u.deg).value)
    fix_tab.replace_column(
        "leakage",
        Column(
            leakage,
            name="leakage",
        ),
    )

    goodI = ~fix_tab["stokesI_fit_flag"] & ~fix_tab["channel_flag"]
    goodL = goodI & ~fix_tab["leakage_flag"] & (fix_tab["snr_polint"] > 5)
    goodRM = goodL & ~fix_tab["snr_flag"]
    good_fix_tab = fix_tab[goodRM]
    fix_flag_tab = compute_local_rm_flag(good_cat=good_fix_tab, big_cat=fix_tab)

    _, ext = os.path.splitext(cat)
    outfile = cat.replace(ext, f".corrected{ext}")

    outfit = cat.replace(ext, ".corrected.leakage.pkl")
    with open(outfit, "wb") as f:
        pickle.dump(fit, f)
        logger.info(f"Wrote leakage fit to {outfit}")

    logger.info(f"Writing corrected catalogue to {outfile}")
    if ext == ".xml" or ext == ".vot":
        write_votable(fix_flag_tab, outfile)
    else:
        tab.write(outfile, overwrite=True)
    logger.info(f"{outfile} written to disk")
    logger.info("Done!")


def cli():
    import argparse

    parser = argparse.ArgumentParser(description="Fix DR1 catalogs")
    parser.add_argument("catalogue", type=str, help="Input catalog")
    parser.add_argument(
        "survey",
        type=str,
        help="Survey directory",
    )
    parser.add_argument(
        "--epoch",
        type=int,
        default=0,
        help="Epoch to read field data from",
    )
    parser.add_argument("--debug", action="store_true", help="Print debug messages")
    args = parser.parse_args()

    logger.setLevel(logging.INFO)

    if args.debug:
        logger.setLevel(logging.DEBUG)
    main(
        cat=args.catalogue,
        survey_dir=Path(args.survey),
        epoch=args.epoch,
    )


if __name__ == "__main__":
    cli()
