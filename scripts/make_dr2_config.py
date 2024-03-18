#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Make a DR2 config file"""

import logging
import os
from pathlib import Path

import pandas as pd

from arrakis.logger import logger

logger.setLevel(logging.INFO)
RACSUSER = "anonymous"
RACSPASS = "racs-db-letmein"
RACSHOST = "146.118.68.63"
RACSPORT = "5433"
RACSDB = "epoch_7"


def get_field_data(sbid: int) -> pd.Series:
    """Get field data from the RACS database

    Args:
        sbid (int): SBID

    Returns:
        pd.Series: Field data row
    """
    table = "field_data"
    df = pd.read_sql(
        f"SELECT * from {table}",
        f"postgresql://{RACSUSER}:{RACSPASS}@{RACSHOST}:{RACSPORT}/{RACSDB}",
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
    # From Uncle Timmy:
    # SBID 38307-38528 : associated holography is SBID 38585
    # SBID 38545-39385 : associated holography is SBID 38709
    # SBID 39400-40878 : associated holography is SBID 39549
    # SBID 40989-41829 : associated holography is SBID 41055
    holo_dir = Path("/scratch3/projects/spiceracs/RACS_Low2_Holography")
    if 38307 <= sbid <= 38528:
        holo_sbid = 38585
    elif 38545 <= sbid <= 39385:
        holo_sbid = 38709
    elif 39400 <= sbid <= 40878:
        holo_sbid = 39549
    elif 40989 <= sbid <= 41829:
        holo_sbid = 41055
    else:
        raise ValueError(f"SBID {sbid} not in range")
    holo_file = holo_dir / f"akpb.iquv.square_6x6.63.887MHz.SB{holo_sbid}.cube.fits"
    assert holo_file.exists(), f"{holo_file} does not exist"
    return holo_file


def main(
    sbid: int,
    processing_dir: Path,
):
    """Main script"""
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
        ms_glob_pattern="*beam[0-36]*.ms",
        imager_dask_config="/scratch3/projects/spiceracs/arrakis/arrakis/configs/petrichor.yaml",
        mgain=0.7,
        force_mask_rounds=8,
        nmiter=15,
        niter=500_000,
        local_rms=True,
        auto_mask=3.5,
        local_rms_window=60,
        auto_threshold=1,
        size=6144,
        scale=2.5,
        robust=-0.5,
        pols="IQU",
        gridder="wgridder",
        minuv=200,
        local_wsclean="/datasets/work/sa-mhongoose/work/containers/wsclean_force_mask.sif",
        multiscale=False,
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
        fitRMSF=True,
        fit_function="log",
        gain=0.1,
        holofile=holo_file.as_posix(),
        maxIter=10000,
        nSamples=100.0,
        noStokesI=False,
        not_RMSF=False,
        outfile=f"{field_data.FIELD_NAME}_SB{sbid}_polcat.fits",
        pad=5.0,
        polyOrd=-2,
        rm_verbose=False,
        savePlots=True,
        showPlots=False,
        validate=False,
        verbose=True,
        weightType="variance",
        yanda_image="/datasets/work/sa-mhongoose/work/containers/askapsoft_1.15.0-openmpi4.sif",
        ionex_server="file:///datasets/work/sa-mhongoose/work/data/IONEX/ftp.aiub.unibe.ch",
        ionex_formatter="ftp.aiub.unibe.ch",
        ionex_prefix="codg",
    )

    config_file = processing_dir / f"{sbid}_rm.cfg"

    with open(config_file, "w") as f:
        for key, value in config_base.items():
            f.write(f"{key} = {value}\n")

    # Now make a run script
    script_file = processing_dir / f"{sbid}_rm_run.sh"
    script_string = f"""#!/bin/bash
#SBATCH --job-name=spice_master
#SBATCH --export=NONE
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=36GB
#SBATCH --time=1-00:00:00
#SBATCH -A OD-217087
#SBATCH -o {sbid}_rm_%j.log
#SBATCH -e {sbid}_rm_%j.log
#SBATCH --qos=express

# I trust nothing
export OMP_NUM_THREADS=1

export APIURL=http://stokes.it.csiro.au:4200/api
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
	--config {config_file.absolute.as_posix()}  \
	{processing_dir/sbid} \
	{processing_dir} \
	{field_data.FIELD_NAME} \
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
        "-p", "--procdir", type=Path, help="Processing directory", default=Path(".")
    )

    args = parser.parse_args()

    _ = main(args.sbid, args.procdir)


if __name__ == "__main__":
    cli()
