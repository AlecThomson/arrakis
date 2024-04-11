#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Make a DR2 config file"""

import logging
import os
from pathlib import Path

import pandas as pd
import yaml

from arrakis.logger import logger

logger.setLevel(logging.INFO)
# Public RACS database credentials
RACSUSER = "anonymous"
RACSHOST = "146.118.68.63"
RACSPORT = "5433"


def get_field_data(sbid: int) -> pd.Series:
    """Get field data from the RACS database

    Args:
        sbid (int): SBID

    Returns:
        pd.Series: Field data row
    """
    table = "field_data"

    if 55538 <= sbid <= 59072:
        racsdb = "epoch_9"
    elif 38307 <= sbid <= 41829:
        racsdb = "epoch_7"

    df = pd.read_sql(
        f"SELECT * from {table}",
        f"postgresql://{RACSUSER}:{os.environ['PGPASSWORD']}@{RACSHOST}:{RACSPORT}/{racsdb}",
    )
    df.set_index("SBID", inplace=True)
    return df.loc[sbid]


def get_holography_path(sbid: int) -> Path:
    """Get the path to the holography file for a given SBID

    Args:
        sbid (int): SBID

    Raises:
        ValueError: If SBID is not in the expected range

    Returns:
        Path: Path to the holography file
    """

    # RACS-low3 min SBID 55538
    # SBID 55538-59072 : associated holography is SBID 55219

    # RACS-low2 38307 - 41829
    # From Uncle Timmy:
    # SBID 38307-38528 : associated holography is SBID 38585
    # SBID 38545-39385 : associated holography is SBID 38709
    # SBID 39400-40878 : associated holography is SBID 39549
    # SBID 40989-41829 : associated holography is SBID 41055

    if 55538 <= sbid <= 59072:
        holo_dir = Path("/scratch3/projects/spiceracs/RACS_Low3_Holography")
    elif 38307 <= sbid <= 41829:
        holo_dir = Path("/scratch3/projects/spiceracs/RACS_Low2_Holography")
    else:
        raise ValueError(f"SBID {sbid} not in range")

    if 38307 <= sbid <= 38528:
        holo_sbid = 38585
    elif 38545 <= sbid <= 39385:
        holo_sbid = 38709
    elif 39400 <= sbid <= 40878:
        holo_sbid = 39549
    elif 40989 <= sbid <= 41829:
        holo_sbid = 41055
    elif 55538 <= sbid <= 59072:
        holo_sbid = 55219
    else:
        raise ValueError(f"SBID {sbid} not in range")

    if 38307 <= sbid <= 41829:
        holo_file = holo_dir / f"akpb.iquv.square_6x6.63.887MHz.SB{holo_sbid}.cube.fits"

    elif 55538 <= sbid <= 59072:
        holo_file = (
            holo_dir / f"akpb.iquv.closepack36.54.943MHz.SB{holo_sbid}.cube.fits"
        )

    assert holo_file.exists(), f"{holo_file} does not exist"
    return holo_file


def main(
    sbid: int,
    sbid_dir: Path,
    processing_dir: Path,
):
    """Main script"""

    if not processing_dir.exists():
        processing_dir.mkdir(parents=True, exist_ok=True)

    field_data = get_field_data(sbid)
    holo_file = get_holography_path(sbid)

    # Set nchan depending on the Galactic latitude
    nchan = 36 if abs(field_data.GAL_LAT) > 10 else 72

    config_base = dict(
        host="stokes.it.csiro.au",
        username="admin",
        password=os.environ["SPICE_PASSWORD"],
        imager_only=False,
        epoch=7,
        ms_glob_pattern=f"'SB{sbid}.{field_data.FIELD_NAME}.beam*.round4.ms'",
        imager_dask_config="/scratch3/projects/spiceracs/arrakis/arrakis/configs/petrichor.yaml",
        temp_dir_images="$LOCALDIR",
        temp_dir_wsclean="$MEMDIR",
        mgain=0.7,
        force_mask_rounds=8,
        nmiter=15,
        niter=500_000,
        local_rms=True,
        auto_mask=4,
        local_rms_window=60,
        auto_threshold=1,
        size=6144,
        scale=2.5,
        robust=-0.5,
        pols="IQU",
        gridder="wgridder",
        minuv=200,
        local_wsclean="/datasets/work/sa-mhongoose/work/containers/wsclean_force_mask.sif",
        multiscale=True,
        multiscale_scale_bias=0.6,
        multiscale_scales="0,2,4,8,16,32,64,128",
        purge=True,
        absmem=100,
        nchan=nchan,
        skip_fix_ms=True,
        data_column="CORRECTED_DATA",
        skip_imager=False,
        skip_cutout=False,
        skip_linmos=False,
        skip_cleanup=False,
        skip_frion=False,
        skip_rmsynth=False,
        skip_rmclean=False,
        skip_cat=False,
        cutoff=-8.0,
        window=-5.0,
        dask_config="/scratch3/projects/spiceracs/arrakis/arrakis/configs/rm_petrichor.yaml",
        database=True,
        debug=False,
        dimension="1d",
        dryrun=False,
        fit_rmsf=True,
        fit_function="log",
        gain=0.1,
        holofile=holo_file.as_posix(),
        max_iter=10000,
        n_samples=100.0,
        no_stokes_i=False,
        not_rmsf=False,
        write=f"{field_data.FIELD_NAME}_SB{sbid}_polcat.fits",
        pad=7.0,
        poly_ord=-2,
        rm_verbose=False,
        save_plots=True,
        show_plots=False,
        validate=False,
        weight_type="variance",
        yanda_image="/datasets/work/sa-mhongoose/work/containers/askapsoft_1.15.0-openmpi4.sif",
        ionex_server="file:///datasets/work/sa-mhongoose/work/data/IONEX/gdc.cddis.eosdis.nasa.gov",
        ionex_formatter="cddis.gsfc.nasa.gov",
        ionex_prefix="casg",
    )

    config_file = processing_dir / f"{sbid}_rm.yaml"

    with open(config_file, "w") as f:
        yaml.safe_dump(config_base, f)

    # Now make a run script
    script_file = processing_dir / f"{sbid}_rm_run.sh"
    script_string = rf"""#!/bin/bash
#SBATCH --job-name=spice_master
#SBATCH --export=NONE
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=36GB
#SBATCH --time=1-00:00:00
#SBATCH -A OD-217087
#SBATCH -o {(processing_dir/str(sbid)).absolute().as_posix()}_rm_%j.log
#SBATCH -e {(processing_dir/str(sbid)).absolute().as_posix()}_rm_%j.log
#SBATCH --qos=express

# I trust nothing
export OMP_NUM_THREADS=1

export APIURL=http://jones.it.csiro.au:4200/api
export PREFECT_API_URL="${{APIURL}}"
export WORKDIR=$(pwd)
export PREFECT_HOME="${{WORKDIR}}/prefect"
export PREFECT_LOGGING_EXTRA_LOGGERS="arrakis"

echo "Sourcing home"
source /home/$(whoami)/.bashrc

module load singularity

echo "Activating conda arrakis environment"
conda activate arrakis310

echo "About to run spice_process"
spice_process \
	{processing_dir.absolute().as_posix()} \
	{field_data.FIELD_NAME} \
	{sbid_dir.absolute().as_posix()} \
    --sbid {sbid} \
	--config {config_file.absolute().as_posix()}  \
"""
    with open(script_file, "w") as f:
        f.write(script_string)

    logger.info(f"Wrote {config_file} and {script_file}")

    return config_file, script_file


def cli():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("sbid", type=int, help="SBID")
    parser.add_argument(
        "-s", "--sbiddir", type=Path, help="Processing directory", default=Path(".")
    )
    parser.add_argument(
        "-p", "--procdir", type=Path, help="Processing directory", default=Path(".")
    )

    args = parser.parse_args()

    _ = main(sbid=args.sbid, sbid_dir=args.sbiddir, processing_dir=args.procdir)


if __name__ == "__main__":
    cli()
