#!/usr/bin/env python3
from astropy.table import Table, Row
import argparse
import os


def main(name: str, cal=False, science=False):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    basedir = f"{scriptdir}/../askap_surveys/racs/db/epoch_0"
    tab = Table.read(f'{basedir}/field_data.csv')
    tab.add_index('FIELD_NAME')
    sel_tab = tab[tab['SELECT'] == 1]
    sub_tab = Table(sel_tab.loc['FIELD_NAME', f"RACS_{name}"])
    space = '       '
    if cal:
        print(int(sub_tab['CAL_SBID']))

    if science:
        print(int(sub_tab['SBID']))
    if not cal and not science:
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
    parser.add_argument(
        "--cal",
        action='store_true',
        help="Return CAL SBID only"
    )
    parser.add_argument(
        "--science",
        action='store_true',
        help="Return Science SBID only"
    )
    args = parser.parse_args()
    main(
        args.field,
        cal=args.cal,
        science=args.science
    )
