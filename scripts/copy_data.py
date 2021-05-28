#!/usr/bin/env python3
from astropy.table import Table
import argparse
import os
from glob import glob
from spiceracs.utils import try_mkdir
from shutil import copyfile, SameFileError


def yes_or_no(question):
    while "Please answer 'y' or 'n'":
        reply = str(input(question+' (y/n): ')).lower().strip()
        if reply[:1] == 'y':
            return True
        if reply[:1] == 'n':
            return False


def rsync(src, tgt):
    os.system(f"rsync -rPvh {src} {tgt}")


def prsync(wild_src: str, tgt: str, ncores: int):
    os.system(
        f"ls -d {wild_src} | xargs -n1 -P{ncores} -I% rsync -rPvh % {tgt}")


def main(name: str,
         sbid: int,
         racs_area: str,
         spice_area: str,
         ncores=1,
         clean=False,
         force=False
         ):
    tab = Table.read(
        '/group/askap/athomson/repos/spiceracs/askap_surveys/racs/db/epoch_0/field_data.csv'
    )
    tab.add_index('FIELD_NAME')
    tab.add_index('CAL_SBID')
    row = Table(
        tab.loc['FIELD_NAME', f"RACS_{name}"]
    ).loc['CAL_SBID', sbid]['INDEX']
    sb_dir = os.path.abspath(f"{spice_area}/{sbid}")
    field_dir = os.path.abspath(f"{sb_dir}/RACS_test4_1.05_{name}")
    bpcal = os.path.abspath(f"{sb_dir}/BPCAL")
    check = os.path.abspath(f"{field_dir}/Checkfiles")

    for idir in [sb_dir, field_dir, bpcal, check]:
        try_mkdir(idir)

    prsync(
        f"{racs_area}/{sbid}/BPCAL/calparameters_1934_bp_*.tab",
        f"{bpcal}/",
        ncores
    )
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
        ncores
    )

    rsync(
        f"{racs_area}/{sbid}/RACS_test4_1.05_{name}/Checkfiles/",
        f"{check}/"
    )

    # Fix for multiple fields
    for f in sorted(glob(f'{check}/*')):
        abspath = os.path.abspath(f)
        idx = abspath.find('_F')
        f_no = abspath[idx+1:idx+4]
        newpath = abspath.replace(f_no, 'F00')
        try:
            copyfile(abspath, newpath)
        except SameFileError:
            pass
        print(os.path.basename(newpath))

    if clean:
        if force:
            yes = True
        else:
            yes = yes_or_no(
                f"This will delete the CONTCUBE checkfiles in {check}. Are you sure?")
        if yes:
            files = glob(f"{check}/CONTCUBE*")
            for f in files:
                os.remove(f)


if __name__ == "__main__":
    descStr = f"""
    Copy data from RACS area to SPICE area'
    """
    parser = argparse.ArgumentParser(
        description=descStr,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "field",
        metavar="field",
        type=str,
        help="RACS field to find e.g. 2132-50A"
    )

    parser.add_argument(
        "cal_sbid",
        metavar="cal_sbid",
        type=int,
        help="Calibrator SBID for field",
    )
    parser.add_argument(
        "--clean",
        action='store_true',
        help="Cleanup Checkfiles",
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
        default=os.path.abspath('/askapbuffer/payne/mcc381/RACS'),
        help="RACS area",
    )
    parser.add_argument(
        "--spice",
        type=str,
        default=os.path.abspath(
            '/group/askap/athomson/projects/spiceracs/spica'),
        help="SPICE area",
    )

    args = parser.parse_args()
    main(args.field, args.cal_sbid, args.ncores, clean=args.clean)
