#!/usr/bin/env python3
import argparse
import os
from glob import glob
from pathlib import Path
from shutil import SameFileError, copyfile

from astropy.table import Table

from arrakis.logger import logger
from arrakis.utils import try_mkdir


def yes_or_no(question):
    while "Please answer 'y' or 'n'":
        reply = str(input(question + " (y/n): ")).lower().strip()
        if reply[:1] == "y":
            return True
        if reply[:1] == "n":
            return False


def rsync(src, tgt):
    os.system(f"rsync -rPvh {src} {tgt}")


def prsync(wild_src: str, tgt: str, ncores: int):
    os.system(f"ls -d {wild_src} | xargs -n 1 -P {ncores} -I% rsync -rvh % {tgt}")


def main(
    name: str,
    sbid: int,
    racs_area: str,
    spice_area: str,
    survey_dir: Path,
    epoch: int = 0,
    ncores=1,
    clean=False,
    force=False,
):
    field_path = survey_dir / "db" / f"epoch_{epoch}" / "field_data.csv"
    tab = Table.read(field_path)
    tab.add_index("FIELD_NAME")
    tab.add_index("CAL_SBID")
    row = Table(tab.loc["FIELD_NAME", f"RACS_{name}"]).loc["CAL_SBID", sbid]["INDEX"]
    sb_dir = os.path.abspath(f"{spice_area}/{sbid}")
    field_dir = os.path.abspath(f"{sb_dir}/RACS_test4_1.05_{name}")
    bpcal = os.path.abspath(f"{sb_dir}/BPCAL")
    check = os.path.abspath(f"{field_dir}/Checkfiles")
    if clean:
        if force:
            yes = True
        else:
            yes = yes_or_no(
                f"This will delete the CONTCUBE checkfiles in {check}. Are you sure?"
            )
    for idir in [sb_dir, field_dir, bpcal, check]:
        try_mkdir(idir)

    prsync(f"{racs_area}/{sbid}/BPCAL/calparameters_1934_bp_*.tab", f"{bpcal}/", ncores)
    # Needed until pipeline update
    # prsync(
    #     f"{racs_area}/{sbid}/RACS_test4_1.05_{name}/*_averaged_cal.ms",
    #     f"{field_dir}/",
    #     ncores
    # )
    # prsync(
    #     f"{racs_area}/{sbid}/RACS_test4_1.05_{name}/*.ms.flagSummary",
    #     f"{field_dir}/",
    #     ncores
    # )
    prsync(
        f"{racs_area}/{sbid}/RACS_test4_1.05_{name}/*_averaged_cal.leakage.ms",
        f"{field_dir}/",
        ncores,
    )

    rsync(f"{racs_area}/{sbid}/RACS_test4_1.05_{name}/Checkfiles/", f"{check}/")

    # Fix for multiple fields
    for f in sorted(glob(f"{check}/*")):
        abspath = os.path.abspath(f)
        idx = abspath.find("_F")
        f_no = abspath[idx + 1 : idx + 4]
        newpath = abspath.replace(f_no, "F00")
        try:
            copyfile(abspath, newpath)
        except SameFileError:
            pass
        logger.debug(os.path.basename(newpath))

    if clean:
        if yes:
            files = glob(f"{check}/CONTCUBE*")
            for f in files:
                os.remove(f)


def cli():
    descStr = f"""
    Copy data from RACS area to SPICE area'
    """
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field to find e.g. 2132-50A"
    )

    parser.add_argument(
        "cal_sbid",
        metavar="cal_sbid",
        type=int,
        help="Calibrator SBID for field",
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
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Cleanup Checkfiles",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force cleanup of Checkfiles",
    )
    parser.add_argument(
        "--ncores",
        type=int,
        default=1,
        help="Ncores for parallel rsync",
    )

    parser.add_argument(
        "--RACS",
        type=str,
        default=os.path.abspath("/askapbuffer/payne/mcc381/RACS"),
        help="RACS area",
    )
    parser.add_argument(
        "--spice",
        type=str,
        default=os.path.abspath("/scratch/ja3/athomson/spica"),
        help="SPICE area",
    )

    args = parser.parse_args()
    main(
        name=args.field,
        sbid=args.cal_sbid,
        racs_area=args.RACS,
        spice_area=args.spice,
        survey_dir=Path(args.survey),
        epoch=args.epoch,
        ncores=args.ncores,
        clean=args.clean,
        force=args.force,
    )


if __name__ == "__main__":
    cli()
