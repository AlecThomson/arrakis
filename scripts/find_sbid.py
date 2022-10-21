#!/usr/bin/env python3
import argparse
import os

from astropy.table import Row, Table

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


def main(name: str, cal=False, science=False, weight=False):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    basedir = f"{scriptdir}/../askap_surveys/racs/db/epoch_0"
    tab = Table.read(f"{basedir}/field_data.csv")
    tab.add_index("FIELD_NAME")
    sel_tab = tab[tab["SELECT"] == 1]
    sub_tab = Table(sel_tab.loc["FIELD_NAME", f"RACS_{name}"])
    space = "       "
    if cal:
        print(int(sub_tab["CAL_SBID"]))

    if science:
        print(int(sub_tab["SBID"]))
    if weight:
        sbid = int(sub_tab["SBID"])
        print(int(sorted_weights[sorted_sbids.index(sbid)]))
    if not cal and not science and not weight:
        print(f"DB info for RACS_{name}:\n")
        for i, row in enumerate(sub_tab):
            print(f"{space}CAL SBID {i+1}: {row['CAL_SBID']}")
            print(f"{space}Science SBID {i+1}: {row['SBID']}\n")


def cli():
    descStr = """
    Find SBID(s) in RACS-DB for given field
    """
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field to find e.g. 2132-50A"
    )
    parser.add_argument("--cal", action="store_true", help="Return CAL SBID only")
    parser.add_argument(
        "--science", action="store_true", help="Return Science SBID only"
    )
    parser.add_argument("--weight", action="store_true", help="Return weight SBID only")
    args = parser.parse_args()
    main(
        args.field,
        cal=args.cal,
        science=args.science,
        weight=args.weight,
    )


if __name__ == "__main__":
    cli()
