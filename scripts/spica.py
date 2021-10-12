import os
import subprocess as sb
import shlex
from glob import glob
from pprint import pprint
from astropy.table import Table
import numpy as np
import copy_data
from spiceracs.utils import try_mkdir

racs_area = os.path.abspath('/askapbuffer/payne/mcc381/RACS')
# spice_area = os.path.abspath('/group/askap/athomson/projects/spiceracs/spica')
spice_area = os.path.abspath('/scratch/ja3/athomson/spica')
group_area = os.path.abspath('/group/ja3/athomson/spica/')

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

def mslist(cal_sb, name):
    # os.system('module unload askapsoft')
    # os.system('module load askapsoft')
    try:
        ms = glob(f"{racs_area}/{cal_sb}/RACS_test4_1.05_{name}/*beam00_*.ms")[0]
        # print('ms',ms)
    except:
        raise Exception(f"Can't find '{racs_area}/{cal_sb}/RACS_test4_1.05_{name}/*beam00_*.ms'")

    mslist_out = sb.run(shlex.split(f"mslist --full {ms}"), capture_output=True, check=False)
    if mslist_out.returncode > 0:
        print(mslist_out.stderr.decode('utf-8'))
        print(mslist_out.stdout.decode('utf-8'))
        mslist_out.check_returncode()
    date_out = sb.run(shlex.split('date +%Y-%m-%d-%H%M%S'), capture_output=True, check=True)

    out = mslist_out.stderr.decode() + f"METADATA_IS_GOOD {date_out.stdout.decode()}"
    # print(out)
    return out

def main(copy=False, force=False, cal=False, mslist_dir=None, cube_image=False):
    tab = Table.read(f'{basedir}/field_data.csv')
    tab.add_index('FIELD_NAME')

    cols = [
        'Field name',
        'CAL SBID',
        'SBID',
        'Leakage cal',
        'Row index',
        'Cube imaging'
    ]
    spica_tab = Table(names=cols, dtype=[str, int, int, bool, int, bool])
    for name in SPICA:
        sub_tab = Table(tab.loc['FIELD_NAME', f"RACS_{name}"])
        idx = np.argmax(sub_tab['CAL_SBID'])
        index = sub_tab['INDEX'][idx]
        cal_sbid = sub_tab['CAL_SBID'][idx]
        sbid = sub_tab['SBID'][idx]
        cal_files = glob(
            f"{racs_area}/{cal_sbid}/RACS_test4_1.05_{name}/*averaged_cal.leakage.ms")
        leak = not len(cal_files) == 0
        cubes = []
        for stoke in ['i', 'q', 'u']:
            cubes.extend(glob(
                f"{spice_area}/{cal_sbid}/RACS_test4_1.05_{name}/image.restored.{stoke}*.conv.fits"))
            cubes.extend(glob(
            f"{spice_area}/{cal_sbid}/RACS_test4_1.05_{name}/weights.{stoke}*.txt"))
        image = len(cubes) == 216
        spica_tab.add_row([f"RACS_{name}", cal_sbid, sbid, leak, index, image])

    spica_tab.sort('SBID')
    spica_tab.pprint_all()
    if cal:
        print('The following row indcies are ready to image:')
        sub_tab = spica_tab[spica_tab['Leakage cal']]
        indxs = []
        for row in sub_tab:
            indxs.append(row['Row index'])
        indxs = np.array(indxs)
        print(' '.join(indxs.astype(str)))

    if mslist_dir is not None:
        mslist_dir = os.path.abspath(mslist_dir)
        # print('mslist_dir',mslist_dir)
        for row in spica_tab:
            try:
                out = mslist(
                    name=row['Field name'].replace('RACS_', ''),
                    cal_sb=row['CAL SBID']
                )
                # print('out',out)
                sbdir = f"{mslist_dir}/{row['SBID']}"
                # print('sbdir',sbdir)
                try_mkdir(sbdir, verbose=False)
                outdir = f"{sbdir}/metadata"
                try_mkdir(outdir, verbose=False)
                outfile = f"{outdir}/mslist-20_{row['Field name']}.txt"
                with open(outfile, 'w') as f:
                    f.write(out)
            except Exception as e:
                print(e)
                continue

    if copy:
        for row in spica_tab:
            if row['Leakage cal']:
                if row['Cube imaging']:
                    print(f"Cube imaging done for {row['Field name']}. Skipping...")
                    continue
                else:
                    copy_data.main(
                        name=row['Field name'].replace('RACS_', ''),
                        sbid=row['CAL SBID'],
                        racs_area=racs_area,
                        spice_area=spice_area,
                        ncores=10,
                        clean=True,
                        force=force
                    )

    if cube_image:
        cmds = [f"cd {spice_area}", "conda activate aces"]
        for row in spica_tab:
            cmd = f"start_pipeline.py -e 0 -p /group/askap/athomson/projects/spiceracs/spica/racs_pipeline_cube.parset -o -m /group/askap/athomson/projects/spiceracs/spica/modules.txt -t /group/askap/athomson/projects/spiceracs/MSlists/{row['SBID']}/metadata/ -i {row['SBID']} -c {row['CAL SBID']} -f {row['Field name'].replace('RACS_', 'RACS_test4_1.05_')} -a ja3 -s"
            cmds.append(cmd)
        print(f"Written imaging commands to '{cube_image}'")
        with open(cube_image, 'w') as f:
            f.write("\n".join(cmds))
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
        "--copy_cutouts",
        action="store_true",
        help="Copy cutouts back to /group [False]."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force cleanup of Checkfiles."
    )
    parser.add_argument(
        "--cal",
        action="store_true",
        help="Print calibrated field row indices."
    )
    parser.add_argument(
        "--mslist_dir",
        type=str,
        default=None,
        help="Dir to store mslist files."
    )

    parser.add_argument(
        "--cube_image",
        type=str,
        default=False,
        help="File to write image cmds to."
    )

    args = parser.parse_args()
    _ = main(
        copy=args.copy,
        force=args.force,
        cal=args.cal,
        mslist_dir=args.mslist_dir,
        cube_image=args.cube_image
    )
