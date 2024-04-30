import os
import shlex
import subprocess as sb
from glob import glob
from pathlib import Path

import copy_data
import numpy as np
from astropy.table import Table

from arrakis.logger import logger, logging
from arrakis.utils.io import try_mkdir

logger.setLevel(logging.INFO)

racs_area = os.path.abspath("/askapbuffer/payne/mcc381/RACS")
# spice_area = os.path.abspath('/group/askap/athomson/projects/arrakis/spica')
spice_area = os.path.abspath("/scratch/ja3/athomson/spica")
group_area = os.path.abspath("/group/ja3/athomson/spica/")

SPICA = [
    "1416+00A",
    "1351+00A",
    "1326+00A",
    "1237+00A",
    "1302+00A",
    "1416-06A",
    "1351-06A",
    "1326-06A",
    "1302-06A",
    "1237-06A",
    "1418-12A",
    "1353-12A",
    "1328-12A",
    "1303-12A",
    "1237-12A",
    "1424-18A",
    "1357-18A",
    "1331-18A",
    "1305-18A",
    "1239-18A",
    "1429-25A",
    "1402-25A",
    "1335-25A",
    "1307-25A",
    "1240-25A",
    "1212+00A",
    "1212-06A",
    "1212-12A",
    "1213-18A",
    "1213-25A",
]


def mslist(cal_sb, name):
    try:
        ms = glob(f"{racs_area}/{cal_sb}/RACS_test4_1.05_{name}/*beam00_*.ms")[0]
    except Exception as e:
        logger.error(e)
        raise Exception(
            f"Can't find '{racs_area}/{cal_sb}/RACS_test4_1.05_{name}/*beam00_*.ms'"
        )

    mslist_out = sb.run(
        shlex.split(f"mslist --full {ms}"), capture_output=True, check=False
    )
    if mslist_out.returncode > 0:
        logger.error(mslist_out.stderr.decode("utf-8"))
        logger.error(mslist_out.stdout.decode("utf-8"))
        mslist_out.check_returncode()
    date_out = sb.run(
        shlex.split("date +%Y-%m-%d-%H%M%S"), capture_output=True, check=True
    )

    out = mslist_out.stderr.decode() + f"METADATA_IS_GOOD {date_out.stdout.decode()}"
    return out


def main(
    survey_dir: Path,
    epoch: int = 0,
    copy=False,
    force=False,
    cal=False,
    mslist_dir=None,
    cube_image=False,
):
    field_path = survey_dir / "db" / f"epoch_{epoch}" / "field_data.csv"
    tab = Table.read(field_path)
    tab.add_index("FIELD_NAME")

    cols = [
        "Field name",
        "CAL SBID",
        "SBID",
        "Leakage cal",
        "Row index",
        "Cube imaging",
    ]
    spica_tab = Table(names=cols, dtype=[str, int, int, bool, int, bool])
    for name in SPICA:
        sub_tab = Table(tab.loc["FIELD_NAME", f"RACS_{name}"])
        idx = np.argmax(sub_tab["CAL_SBID"])
        index = sub_tab["INDEX"][idx]
        cal_sbid = sub_tab["CAL_SBID"][idx]
        sbid = sub_tab["SBID"][idx]
        cal_files = glob(
            f"{racs_area}/{cal_sbid}/RACS_test4_1.05_{name}/*averaged_cal.leakage.ms"
        )
        leak = not len(cal_files) == 0
        cubes = []
        for stoke in ["i", "q", "u"]:
            cubes.extend(
                glob(
                    f"{spice_area}/{cal_sbid}/RACS_test4_1.05_{name}/image.restored.{stoke}*.conv.fits"
                )
            )
            cubes.extend(
                glob(
                    f"{spice_area}/{cal_sbid}/RACS_test4_1.05_{name}/weights.{stoke}*.txt"
                )
            )
        image = len(cubes) == 216
        spica_tab.add_row([f"RACS_{name}", cal_sbid, sbid, leak, index, image])

    spica_tab.sort("SBID")
    spica_tab.pprint_all()
    if cal:
        logger.info("The following row indcies are ready to image:")
        sub_tab = spica_tab[spica_tab["Leakage cal"]]
        idx_list = []
        for row in sub_tab:
            idx_list.append(row["Row index"])
        indxs = np.array(idx_list)
        logger.info(" ".join(indxs.astype(str)))

    if mslist_dir is not None:
        mslist_dir = os.path.abspath(mslist_dir)
        for row in spica_tab:
            try:
                out = mslist(
                    name=row["Field name"].replace("RACS_", ""), cal_sb=row["CAL SBID"]
                )
                sbdir = f"{mslist_dir}/{row['SBID']}"
                try_mkdir(sbdir, verbose=False)
                outdir = f"{sbdir}/metadata"
                try_mkdir(outdir, verbose=False)
                outfile = f"{outdir}/mslist-20_{row['Field name']}.txt"
                with open(outfile, "w") as f:
                    f.write(out)
            except Exception as e:
                logger.error(e)
                continue

    if copy:
        for row in spica_tab:
            if row["Leakage cal"]:
                if row["Cube imaging"]:
                    logger.info(
                        f"Cube imaging done for {row['Field name']}. Skipping..."
                    )
                    continue
                else:
                    copy_data.main(
                        name=row["Field name"].replace("RACS_", ""),
                        sbid=row["CAL SBID"],
                        racs_area=racs_area,
                        spice_area=spice_area,
                        survey_dir=survey_dir,
                        epoch=epoch,
                        ncores=10,
                        clean=True,
                        force=force,
                    )

    if cube_image:
        cmds = [f"cd {spice_area}", "conda activate aces"]
        for row in spica_tab:
            cmd = f"start_pipeline.py -e 0 -p /group/askap/athomson/projects/arrakis/spica/racs_pipeline_cube.parset -o -m /group/askap/athomson/projects/arrakis/spica/modules.txt -t /group/askap/athomson/projects/arrakis/MSlists/{row['SBID']}/metadata/ -i {row['SBID']} -c {row['CAL SBID']} -f {row['Field name'].replace('RACS_', 'RACS_test4_1.05_')} -a ja3 -s"
            cmds.append(cmd)
        logger.info(f"Written imaging commands to '{cube_image}'")
        with open(cube_image, "w") as f:
            f.write("\n".join(cmds))
    return spica_tab


def cli():
    import argparse

    descStr = """
    Helper for Spica pilot
    """
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.ArgumentDefaultsHelpFormatter
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
        "--copy",
        action="store_true",
        help="Copy calibrated data from racs's area [False].",
    )
    parser.add_argument(
        "--copy_cutouts",
        action="store_true",
        help="Copy cutouts back to /group [False].",
    )
    parser.add_argument(
        "--force", action="store_true", help="Force cleanup of Checkfiles."
    )
    parser.add_argument(
        "--cal", action="store_true", help="Print calibrated field row indices."
    )
    parser.add_argument(
        "--mslist_dir", type=str, default=None, help="Dir to store mslist files."
    )

    parser.add_argument(
        "--cube_image", type=str, default=False, help="File to write image cmds to."
    )

    args = parser.parse_args()
    _ = main(
        survey_dir=Path(args.survey),
        epoch=args.epoch,
        copy=args.copy,
        force=args.force,
        cal=args.cal,
        mslist_dir=args.mslist_dir,
        cube_image=args.cube_image,
    )


if __name__ == "__main__":
    cli()
