#!/usr/bin/env python3
"""Post process DR1 catalog"""
import logging as log
import os
import pickle

import astropy.units as u
import numpy as np
import pkg_resources
from astropy.coordinates import SkyCoord
from astropy.table import Column, Table
from astropy.time import Time
from IPython import embed
from rmtable import RMTable
from spica import SPICA

from spiceracs.makecat import (
    compute_local_rm_flag,
    get_fit_func,
    is_leakage,
    write_votable,
)


def fix_fields(tab: Table) -> Table:
    # Get field data, and index by field/tile ID
    survey_dir = pkg_resources.resource_filename("spiceracs", "askap_surveys")
    basedir = os.path.join(survey_dir, "racs", "db", "epoch_0")
    field = Table.read(os.path.join(basedir, "field_data.csv"))
    field = field[field["SELECT"] == 1]
    field.add_index("FIELD_NAME")
    tab.add_index("tile_id")

    # Compare the fields we have to those we want
    fields_in_cat = list(set(tab["tile_id"]))
    fields_in_spica = [f"RACS_{name}" for name in SPICA]
    log.debug(f"Fields in catalogue: {fields_in_cat}")
    log.debug(f"Fields in spica: {fields_in_spica}")
    fields_not_in_spica = [f for f in fields_in_cat if f not in fields_in_spica]
    spica_field = field.loc[fields_in_spica]
    spica_field_coords = SkyCoord(
        spica_field["RA_DEG"], spica_field["DEC_DEG"], unit=(u.deg, u.deg), frame="icrs"
    )
    start_times = Time(spica_field["SCAN_START"] * u.second, format="mjd")
    spica_field["start_time"] = start_times
    # These are the sources to update
    sources_to_fix = tab.loc[fields_not_in_spica]
    log.info(f"Found {len(sources_to_fix)} sources to fix")

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
    new_tab["tile_id"] = all_fields

    all_seps = (
        new_tab["separation_tile_centre"].value * new_tab["separation_tile_centre"].unit
    )
    all_seps[idx] = min_seps

    all_sbids = new_tab["sbid"].value
    all_sbids[idx] = spica_field["SBID"][min_idx].value

    all_start_times = new_tab["start_time"].value
    all_start_times[idx] = spica_field["start_time"][min_idx].value

    # Update the columns
    new_tab["separation_tile_centre"] = Column(
        data=all_seps,
    )
    new_tab["beamdist"] = Column(
        data=all_seps,
    )

    new_tab["sbid"] = Column(
        data=all_sbids,
    )

    new_tab["start_time"] = Column(
        data=all_start_times,
    )

    # Fix the units - Why does VOTable do this?? Thanks I hate it
    dumb_units = {
        "Jy.beam-1": u.Jy / u.beam,
        "mJy.beam-1": u.mJy / u.beam,
    }
    for col in new_tab.colnames:
        if str(new_tab[col].unit) in dumb_units.keys():
            new_tab[col].unit = dumb_units[str(new_tab[col].unit)]

    return new_tab


def main(cat: str):
    log.debug(f"Reading {cat}")
    tab = RMTable.read(cat)
    log.debug(f"Fixing {cat}")
    fix_tab = fix_fields(tab)
    fit, fig = get_fit_func(fix_tab, do_plot=True, nbins=16, degree=4)
    fig.savefig("leakage_fit_dr1_fix.pdf")
    leakage_flag = is_leakage(
        fix_tab["fracpol"].value, fix_tab["beamdist"].to(u.deg).value, fit
    )
    fix_tab["leakage_flag"] = leakage_flag
    leakage = fit(fix_tab["separation_tile_centre"].to(u.deg).value)
    fix_tab["leakage"] = leakage

    goodI = ~fix_tab["stokesI_fit_flag"] & ~fix_tab["channel_flag"]
    goodL = goodI & ~fix_tab["leakage_flag"] & (fix_tab["snr_polint"] > 5)
    goodRM = goodL & ~fix_tab["snr_flag"]
    good_fix_tab = fix_tab[goodRM]
    fix_flag_tab = compute_local_rm_flag(good_cat=good_fix_tab, big_cat=fix_tab)

    _, ext = os.path.splitext(cat)
    outfile = cat.replace(ext, f".corrected{ext}")

    outfit = cat.replace(ext, f".corrected.leakage.pkl")
    with open(outfit, "wb") as f:
        pickle.dump(fit, f)
        log.info(f"Wrote leakage fit to {outfit}")

    # outplot = cat.replace(ext, f'.corrected.leakage.pdf')
    # log.info(f"Writing leakage plot to {outplot}")
    # fig.savefig(outplot, dpi=300, bbox_inches='tight')
    log.info(f"Writing corrected catalogue to {outfile}")
    if ext == ".xml" or ext == ".vot":
        write_votable(fix_flag_tab, outfile)
    else:
        tab.write(outfile, overwrite=True)
    log.info(f"{outfile} written to disk")
    log.info("Done!")


def cli():
    import argparse

    parser = argparse.ArgumentParser(description="Fix DR1 catalogs")
    parser.add_argument("catalogue", type=str, help="Input catalog")
    parser.add_argument("--debug", action="store_true", help="Print debug messages")
    args = parser.parse_args()

    log.basicConfig(
        level=log.INFO,
        format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )

    if args.debug:
        log.basicConfig(
            level=log.DEBUG,
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            force=True,
        )

    main(cat=args.catalogue)


if __name__ == "__main__":
    cli()
