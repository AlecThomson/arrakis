from __future__ import annotations

import argparse
import os

import copy_data
import spica
from arrakis.logger import logger, logging
from arrakis.utils.io import try_mkdir

logger.setLevel(logging.INFO)

racs_area = os.path.abspath("/askapbuffer/payne/mcc381/RACS")
# spice_area = os.path.abspath('/group/askap/athomson/projects/arrakis/spica')
spice_area = os.path.abspath("/scratch/ja3/athomson/spica")
group_area = os.path.abspath("/group/ja3/athomson/spica")


def main(field, dry_run=False, ncores=10):
    spica_tab = spica.main()
    spica_tab.add_index("Field name")
    row = spica_tab.loc["Field name", f"{field}"]
    cut_dir = os.path.join(
        spice_area, f"{row['CAL SBID']}", f"RACS_test4_1.05_{field}", "cutouts"
    )
    test_cut = os.path.isdir(cut_dir)
    if not test_cut:
        raise FileNotFoundError(cut_dir)
    else:
        logger.info(f"Copying '{cut_dir}'")

    store_dir = os.path.join(
        group_area, f"{row['CAL SBID']}", f"RACS_test4_1.05_{field}", "cutouts"
    )
    try_mkdir(store_dir)
    test_store = os.path.isdir(store_dir)
    if not test_cut:
        raise FileNotFoundError(test_store)
    else:
        logger.info(f"Storing in '{store_dir}'")

    if not dry_run:
        copy_data.prsync(f"{cut_dir}/*", store_dir, ncores=ncores)


def cli():
    descStr = """
    Copy data from RACS area to SPICE area'
    """
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field to find e.g. 2132-50A"
    )

    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Don't launch rsync",
    )

    parser.add_argument(
        "--ncores",
        type=int,
        default=1,
        help="Ncores for parallel rsync",
    )

    args = parser.parse_args()
    main(
        field=args.field,
        dry_run=args.dry_run,
        ncores=args.ncores,
    )


if __name__ == "__main__":
    cli()
