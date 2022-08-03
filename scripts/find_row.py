#!/usr/bin/env python3
from astropy.table import Table
import argparse
import os


def main(name: str, sbid: int):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    basedir = f"{scriptdir}/../askap_surveys/racs/db/epoch_0"
    tab = Table.read(f"{basedir}/field_data.csv")
    tab.add_index("FIELD_NAME")
    tab.add_index("CAL_SBID")
    row = tab.loc["FIELD_NAME", f"RACS_{name}"].loc["CAL_SBID", sbid]["INDEX"]

    print(f"Row in RACS database is {row}")


if __name__ == "__main__":
    descStr = """
    Find row in RACS-DB for given field
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
    args = parser.parse_args()
    main(args.field, args.cal_sbid)
