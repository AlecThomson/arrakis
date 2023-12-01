#!/usr/bin/env python3
import argparse
from pathlib import Path

from astropy.table import Table

from arrakis.logger import logger, logging

logger.setLevel(logging.INFO)


def main(
    name: str,
    sbid: int,
    survey_dir: Path,
    epoch: int = 0,
):
    field_path = survey_dir / "db" / f"epoch_{epoch}" / "field_data.csv"
    tab = Table.read(field_path)
    tab.add_index("FIELD_NAME")
    tab.add_index("CAL_SBID")
    row = tab.loc["FIELD_NAME", f"RACS_{name}"].loc["CAL_SBID", sbid]["INDEX"]

    logger.info(f"Row in RACS database is {row}")


def cli():
    descStr = """
    Find row in RACS-DB for given field
    """
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.ArgumentDefaultsHelpFormatter
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
    args = parser.parse_args()
    main(
        name=args.field,
        sbid=args.cal_sbid,
        survey_dir=Path(args.survey),
        epoch=args.epoch,
    )


if __name__ == "__main__":
    cli()
