#!/usr/bin/env python3
"""Post process DR1 source catalog"""
import os
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from IPython import embed
from tqdm.auto import tqdm

from arrakis.logger import logger
from arrakis.makecat import fix_blank_units, replace_nans, vot

logger.setLevel("DEBUG")


def add_metadata(vo_table: vot.tree.Table, rmtab: Table, filename: str):
    """Add metadata to VO Table for CASDA

    Args:
        vo_table (vot): VO Table object

    Returns:
        vot: VO Table object with metadata
    """
    # Add metadata
    for col_idx, col_name in enumerate(rmtab.colnames):
        col = vo_table.get_first_table().get_field_by_id(col_name)
        meta_idx = col_idx + 1
        if f"TCOMM{meta_idx}" in rmtab.meta:
            logger.info(f"Adding metadata for {col_name}")
            col.description = rmtab.meta[f"TCOMM{meta_idx}"]
        if f"TUCD{meta_idx}" in rmtab.meta:
            logger.info(f"Adding UCD for {col_name}")
            col.ucd = rmtab.meta[f"TUCD{meta_idx}"]
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


def write_votable(rmtab: Table, outfile: str) -> None:
    # Replace bad column names
    fix_columns = {
        "catalog": "catalog_name",
        "interval": "obs_interval",
    }
    # CASDA needs v1.3
    for col_name, new_name in fix_columns.items():
        if col_name in rmtab.colnames:
            rmtab.rename_column(col_name, new_name)
    # Fix blank units
    rmtab = fix_blank_units(rmtab)
    vo_table = vot.from_table(rmtab)
    vo_table.version = "1.3"
    vo_table = add_metadata(vo_table, rmtab, outfile)
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
    for tile_id in tqdm(tile_ids):
        field_tile_idx = field["FIELD_NAME"] == tile_id
        source_tile_idx = source_cut["Tile_ID"] == tile_id
        tile_coords = field_coords[field_tile_idx]
        seps = source_coords[source_tile_idx].separation(tile_coords)
        source_cut["Separation_Tile_Centre"][source_tile_idx] = seps.deg

    # Fix bad columns
    bad_cols = {
        "min": "min_axis",
        "maj": "maj_axis",
        "e_min": "e_min_axis",
        "e_maj": "e_maj_axis",
    }
    for bad_col, good_col in bad_cols.items():
        if bad_col in source_cut.colnames:
            source_cut.rename_column(bad_col, good_col)


    # Write the output cat
    out_pth = spice_cat_pth.parent / ("RACS_DR1_Sources_" + spice_cat_pth.name)
    logger.info(f"Writing corrected catalogue to {out_pth}")
    # source_cut.write(out_pth, format="votable", overwrite=True)
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
