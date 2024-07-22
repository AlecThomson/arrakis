#!/usr/bin/env python3
import argparse
from pathlib import Path

from arrakis.logger import logger, logging
from astropy.table import Table

logger.setLevel(logging.INFO)

sorted_sbids = [
    8570,
    8574,
    8576,
    8584,
    8589,
    8593,
    8674,
    12420,
    12422,
    12423,
    12425,
    12426,
    12427,
    12428,
    12429,
    12430,
    12431,
    12435,
    12493,
    12494,
    12496,
    12497,
    12500,
    12502,
    13587,
    13588,
    13589,
    13591,
    13595,
    13671,
    13672,
    13673,
    13678,
    13743,
    13746,
    13747,
    13749,
]

# Just the list of fields used in DR1
sorted_weights = [
    8247,
    8247,
    8247,
    8247,
    8247,
    8247,
    8669,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    11371,
    13624,
    13624,
    13624,
    13624,
    13624,
    13624,
    13624,
    13624,
]


def main(
    name: str, survey_dir: Path, epoch: int = 0, cal=False, science=False, weight=False
):
    field_path = survey_dir / "db" / f"epoch_{epoch}" / "field_data.csv"
    tab = Table.read(field_path)
    tab.add_index("FIELD_NAME")
    sel_tab = tab[tab["SELECT"] == 1]
    sub_tab = Table(sel_tab.loc["FIELD_NAME", f"RACS_{name}"])
    space = "       "
    if cal:
        logger.info(int(sub_tab["CAL_SBID"]))

    if science:
        logger.info(int(sub_tab["SBID"]))
    if weight:
        sbid = int(sub_tab["SBID"])
        logger.info(int(sorted_weights[sorted_sbids.index(sbid)]))
    if not cal and not science and not weight:
        logger.info(f"DB info for RACS_{name}:\n")
        for i, row in enumerate(sub_tab):
            logger.info(f"{space}CAL SBID {i+1}: {row['CAL_SBID']}")
            logger.info(f"{space}Science SBID {i+1}: {row['SBID']}\n")


def cli():
    descStr = """
    Find SBID(s) in RACS-DB for given field
    """
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field to find e.g. 2132-50A"
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
    parser.add_argument("--cal", action="store_true", help="Return CAL SBID only")
    parser.add_argument(
        "--science", action="store_true", help="Return Science SBID only"
    )
    parser.add_argument("--weight", action="store_true", help="Return weight SBID only")
    args = parser.parse_args()
    main(
        name=args.field,
        survey_dir=Path(args.survey),
        epoch=args.epoch,
        cal=args.cal,
        science=args.science,
        weight=args.weight,
    )


if __name__ == "__main__":
    cli()
