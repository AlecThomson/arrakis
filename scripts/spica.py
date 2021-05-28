from astropy.table import Table
import os
from glob import glob
import numpy as np
import copy_data

racs_area = os.path.abspath('/askapbuffer/payne/mcc381/RACS')
spice_area = os.path.abspath('/group/askap/athomson/projects/spiceracs/spica')

SPICA = [
    '1416+00A',
    '1351+00A',
    '1326+00A',
    '1237+00A',
    '1302+00A',
    '1416-06A',
    '1351-06A',
    '1326-06A',
    '1302-06A',
    '1237-06A',
    '1418-12A',
    '1353-12A',
    '1328-12A',
    '1303-12A',
    '1237-12A',
    '1424-18A',
    '1357-18A',
    '1331-18A',
    '1305-18A',
    '1239-18A',
    '1429-25A',
    '1402-25A',
    '1335-25A',
    '1307-25A',
    '1240-25A',
    '1212+00A',
    '1212-06A',
    '1212-12A',
    '1213-18A',
    '1213-25A',
]

scriptdir = os.path.dirname(os.path.realpath(__file__))
basedir = f"{scriptdir}/../askap_surveys/racs/db/epoch_0"

def main(copy=False):
    tab = Table.read(f'{basedir}/field_data.csv')
    tab.add_index('FIELD_NAME')

    cols = [
        'Field name',
        'CAL SBID',
        'SBID',
        'Leakage cal'
    ]
    spica_tab = Table(names=cols, dtype=[str, int, int, bool])
    for name in SPICA:
        sub_tab = Table(tab.loc['FIELD_NAME', f"RACS_{name}"])
        idx = np.argmax(sub_tab['CAL_SBID'])
        cal_sbid = sub_tab['CAL_SBID'][idx]
        sbid = sub_tab['SBID'][idx]
        cal_files = glob(f"{racs_area}/{cal_sbid}/RACS_test4_1.05_{name}/*averaged_cal.leakage.ms")
        leak = not len(cal_files) == 0
        spica_tab.add_row([f"RACS_{name}", cal_sbid, sbid, leak])

    print(spica_tab)

    if copy:
        for row in spica_tab:
            if row['Leakage cal']:
                copy_data.main(
                    name=row['Field name'].replace('RACS_', ''),
                    sbid=row['CAL SBID'],
                    racs_area=racs_area,
                    spice_area=spice_area,
                    ncores=10,
                    clean=True,
                    force=True
            )

    return spica_tab

if __name__ == "__main__":
    import argparse

    descStr = """
    Helper for Spica pilot
    """
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--copy",
        action="store_true",
        help="Copy calibrated data from racs's area [False]."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force cleanup of Checkfiles."
    )
    args = parser.parse_args()
    _ = main(args.copy)