#!/usr/bin/env python3
"""Post process DR1 source catalog"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import numpy as np
from arrakis.logger import TqdmToLogger, logger
from arrakis.makecat import fix_blank_units, replace_nans, vot
from astropy.coordinates import SkyCoord
from astropy.table import Table
from tqdm.auto import tqdm

TQDM_OUT = TqdmToLogger(logger, level=logging.INFO)

logger.setLevel("DEBUG")


def add_metadata(vo_table: vot.tree.Table, table: Table, filename: str):
    """Add metadata to VO Table for CASDA

    Args:
        vo_table (vot): VO Table object

    Returns:
        vot: VO Table object with metadata
    """
    # Add metadata
    for col_idx, col_name in enumerate(table.colnames):
        col = table[col_name]
        vocol = vo_table.get_first_table().get_field_by_id(col_name)
        if hasattr(col, "description"):
            logger.info(f"Adding description for {col_name}")
            vocol.description = col.description
        logger.info(f"Adding ucd for {col_name}")
        vocol.ucd = col.meta.get("ucd", "")
    # Add params for CASDA
    if len(vo_table.params) > 0:
        logger.warning(f"{filename} already has params - not adding")
        return vo_table
    _, ext = os.path.splitext(filename)
    cat_name = (
        os.path.basename(filename).replace(ext, "").replace(".", "_").replace("-", "_")
    )
    idx_fields = "Source_ID,Source_Name,Peak_flux,Total_flux_Source,RA,Dec"
    pri_fields = "Source_ID,Source_Name,Peak_flux,Total_flux_Source,RA,Dec"
    params = [
        vot.tree.Param(
            vo_table,
            ID="Catalogue_Name",
            name="Catalogue Name",
            value=cat_name,
            arraysize=str(len(cat_name)),
        ),
        vot.tree.Param(
            vo_table,
            ID="Indexed_Fields",
            name="Indexed Fields",
            value=idx_fields,
            arraysize=str(len(idx_fields)),
        ),
        vot.tree.Param(
            vo_table,
            ID="Principal_Fields",
            name="Principal Fields",
            value=pri_fields,
            arraysize=str(len(pri_fields)),
        ),
    ]
    vo_table.get_first_table().params.extend(params)

    return vo_table


def write_votable(table: Table, outfile: str) -> None:
    # Replace bad column names
    fix_columns = {
        "catalog": "catalog_name",
        "interval": "obs_interval",
    }
    # CASDA needs v1.3
    for col_name, new_name in fix_columns.items():
        if col_name in table.colnames:
            table.rename_column(col_name, new_name)
    # Fix blank units
    table = fix_blank_units(table)
    vo_table = vot.from_table(table)
    vo_table.version = "1.3"
    vo_table = add_metadata(vo_table, table, outfile)
    vot.writeto(vo_table, outfile)
    # Fix NaNs for CASDA
    replace_nans(outfile)


def main(
    source_cat_pth: Path,
    spice_cat_pth: Path,
    survey_dir: Path,
    epoch: int = 0,
):
    logger.info(f"Loading source catalog from {source_cat_pth}")
    source_cat = Table.read(source_cat_pth)
    logger.info(f"Loading SPICE catalog from {spice_cat_pth}")
    spice_cat = Table.read(spice_cat_pth)
    field_path = survey_dir / "db" / f"epoch_{epoch}" / "field_data.csv"
    logger.info(f"Loading field data from {field_path}")
    field = Table.read(field_path)
    field = field[field["SELECT"] == 1]
    field.add_index("FIELD_NAME")

    spice_df = spice_cat.to_pandas()
    spice_df.set_index("source_id", inplace=True)
    spice_df.sort_index(inplace=True)
    source_cat.sort("Source_ID")
    spice_grp = spice_df.groupby("source_id")

    sources = list(set(spice_cat["source_id"]))
    src_idx = np.isin(source_cat["Source_ID"].data.astype(str), sources)
    source_cut = source_cat[src_idx]

    columns_to_fix = {
        "SBID": "sbid",
        "Obs_Start_Time": "start_time",
        "Tile_ID": "tile_id",
    }

    for source_col, spice_col in columns_to_fix.items():
        logger.info(f"Fixing {source_col} with {spice_col}")
        source_cut[source_col] = spice_grp[spice_col].apply(lambda x: x[0])

    # Fix the separations
    logger.info("Fixing separations to field centres")
    field_coords = SkyCoord(field["RA_DEG"], field["DEC_DEG"], unit="deg", frame="icrs")
    source_coords = SkyCoord(
        source_cut["RA"], source_cut["Dec"], unit="deg", frame="icrs"
    )
    tile_ids = np.unique(source_cut["Tile_ID"].data.astype(str))
    for tile_id in tqdm(tile_ids, file=TQDM_OUT):
        field_tile_idx = field["FIELD_NAME"] == tile_id
        source_tile_idx = source_cut["Tile_ID"] == tile_id
        tile_coords = field_coords[field_tile_idx]
        seps = source_coords[source_tile_idx].separation(tile_coords)
        source_cut["Separation_Tile_Centre"][source_tile_idx] = seps.deg

    # Fix bad columns
    bad_cols = {
        "Min": "Min_axis",
        "Maj": "Maj_axis",
        "E_Min": "E_Min_axis",
        "E_Maj": "E_Maj_axis",
    }
    for bad_col, good_col in bad_cols.items():
        if bad_col in source_cut.colnames:
            logger.info(f"Renaming {bad_col} to {good_col}")
            source_cut.rename_column(bad_col, good_col)
        else:
            logger.warning(f"{bad_col} not in source catalog")

    # Write the output cat
    out_pth = spice_cat_pth.parent / ("RACS_DR1_Sources_" + spice_cat_pth.name)
    logger.info(f"Writing corrected catalogue to {out_pth}")
    write_votable(source_cut, out_pth.as_posix())


def cli():
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "source",
        type=str,
        help="RACS Source catalog",
    )
    parser.add_argument(
        "spice",
        type=str,
        help="SPICE-RACS catalog",
    )
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

    args = parser.parse_args()
    main(
        source_cat_pth=Path(args.source),
        spice_cat_pth=Path(args.spice),
        survey_dir=Path(args.survey),
        epoch=args.epoch,
    )


if __name__ == "__main__":
    cli()
