#!/usr/bin/env python3
from astropy.table import Table
import argparse
import os

def main(name: str):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    basedir = f"{scriptdir}/../askap_surveys/racs/db/epoch_0"
    tab = Table.read(f'{basedir}/field_data.csv')
    tab.add_index('FIELD_NAME')
    sub_tab = Table(tab.loc['FIELD_NAME', f"RACS_{name}"])
    space = '       '

    print(f'DB info for RACS_{name}:\n')
    for i, row in enumerate(sub_tab):
        print(f"{space}CAL SBID {i+1}: {row['CAL_SBID']}")
        print(f"{space}Science SBID {i+1}: {row['SBID']}\n")

if __name__ == "__main__":
    descStr = """
    Find SBID(s) in RACS-DB for given field
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
    args = parser.parse_args()
    main(args.field)