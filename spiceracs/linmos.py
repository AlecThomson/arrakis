#!/usr/bin/env /group/askap/athomson/miniconda3/envs/aces/bin/python

import os
import subprocess
import shlex
import ast
import numpy as np
from astropy.table import Table
from glob import glob
from astropy.coordinates import SkyCoord
import astropy.units as u
from tqdm import tqdm, trange
from IPython import embed
from aces.obsplan.config import ACESConfig  # noqa
from astropy.coordinates import SkyCoord
import astropy


def gen_pipline(slurm_dir, fieldname, mongobat, slurmbats, dryrun=True):
    """Generate linmos pipeline

    Args:
        slurm_dir (str): directory for slurm files
        fieldname (str): Name of RACS field
        mongobat (str): Name of mongo batch file
        slurmbats (list): List of linmos batch files
        dryrun (bool, optional): Do a dryrun. Defaults to True.
    """
    pipeline = f"""#! /bin/bash
# first job - no dependencies
RES=$(sbatch {mongobat})
JID=$(echo $RES | cut -d " " -f 4)
"""

    for slurmbat in slurmbats:
        pipeline += f"""
sbatch  --dependency=after:$JID {slurmbat}
    """
    pipefile = f"{slurm_dir}/{fieldname}_linmos_pipeline.sh"
    with open(pipefile, "w") as f:
        f.write(pipeline)
    print(f"Written pipeline file to {pipefile}")
    if not dryrun:
        print(f"Submitting {pipefile}")
        command = f"sh {pipefile}"
        command = shlex.split(command)
        proc = subprocess.run(
            command, encoding="utf-8", check=True
        )


def gen_mongo(fieldname, slurm_dir, logdir):
    """Generate mongodb sbatch

    Args:
        fieldname (str): field name
        slurm_dir (str): directory
        logdir (str): directory
    """
    script = f"""#!/bin/bash -l
#SBATCH --partition=workq
#SBATCH --clusters=galaxy
#SBATCH --exclude=nid00[010,200-202]
#SBATCH --account=askap
# No further constraints applied
# No reservation requested
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --job-name=mongo_linmos
#SBATCH --export=NONE

cd /group/askap/athomson/repos/spiceracs
log={logdir}/{fieldname}_${{SLURM_JOB_ID}}_mongo.log
conda activate py36
numactl --interleave=all mongod --dbpath=database --bind_ip $(hostname -i) --fork --logpath $log
echo $(hostname -i) > {logdir}/mongo_ip.txt
while true; do :; done
wait
"""
    mongobat = f"{slurm_dir}/{fieldname}_mongo.sbatch"
    with open(mongobat, "w") as f:
        f.write(script)
    print(f"Written slurm file to {mongobat}")
    return mongobat


def gen_seps(field):
    """Get beam separations

    Args:
        field (str): File name e.g. 2132-50A

    Returns:
        Table: Separation table
    """
    # get offsets
    aces_cfg = ACESConfig()
    fp_factory = aces_cfg.footprint_factory
    # fp = fp_factory.make_footprint('ak:closepack36', 0.9 * pi / 180., 0.0 * pi / 180.)
    fp = fp_factory.make_footprint(
        'ak:square_6x6', 1.05 * np.pi / 180., 0.0 * np.pi / 180.)

    offsets = np.rad2deg(np.array(fp.offsetsRect))

    master_cat = Table.read('askap_surveys/RACS/admin/epoch_0/field_data.csv')
    master_cat.add_index('FIELD_NAME')
    master_cat = master_cat.loc[f"RACS_test4_1.05_{field}"]
    if type(master_cat) is not astropy.table.row.Row:
        master_cat = master_cat[0]

    # Look for multiple SBIDs - only need one
    cats = glob(
        f'askap_surveys/RACS/admin/epoch_0/beam_inf_*-RACS_test4_1.05_{field}.csv')
    beam_cat = Table.read(cats[0])
    beam_cat.add_index('BEAM_NUM')

    names = ["BEAM", "DELTA_RA", "DELTA_DEC", "BEAM_RA",
             "BEAM_DEC", "FOOTPRINT_RA", "FOOTPRINT_DEC"]

    cols = []
    for beam in range(36):
        beam = int(beam)

        beam_dat = beam_cat.loc[beam]
        beam_coord = SkyCoord(
            beam_dat['RA_DEG']*u.deg, beam_dat['DEC_DEG']*u.deg)
        field_coord = SkyCoord(
            master_cat['RA_DEG']*u.deg, master_cat['DEC_DEG']*u.deg)
        ra_beam = beam_coord.ra.hms
        dec_beam = beam_coord.dec.dms
        ra_field = field_coord.ra.hms
        dec_field = field_coord.dec.dms

        row = [
            beam,
            f"{offsets[beam][0]:0.3f}",
            f"{offsets[beam][1]:0.3f}",
            f"{ra_beam.h:02.0f}:{ra_beam.m:02.0f}:{ra_beam.s:06.3f}",
            f"{dec_beam.d:02.0f}:{abs(dec_beam.m):02.0f}:{abs(dec_beam.s):05.2f}",
            f"{ra_field.h:02.0f}:{ra_field.m:02.0f}:{ra_field.s:06.3f}",
            f"{dec_field.d:02.0f}:{abs(dec_field.m):02.0f}:{abs(dec_field.s):05.2f}",
        ]
        cols += [row]
    tab = Table(data=np.array(cols), names=names)

    return tab


def genparset(field, stoke, datadir, septab, prefix=""):
    """Generate parset for LINMOS

    Args:
        field (str): Name of RACS field
        stoke (str): Stokes parameter
        datadir (str): Directory containing cutouts
        septab (Table): Table of beam seperations
        prefix (str, optional): Search for files with a prefix. Defaults to "".

    Raises:
        Exception: If no files are found in the datadir
    """
    ims = sorted(
        glob(
            f"{datadir}/*.cutout.sm.image.restored.{stoke.lower()}.*.beam*[00-35.fits]"
        )
    )
    if len(ims) == 0:
        print(
            f"{datadir}/*.cutout.sm.image.restored.{stoke.lower()}.*.beam*[00-35.fits]"
        )
        raise Exception(
            'No files found. Have you run imaging? Check your prefix?')
    imlist = "["
    for im in ims:
        imlist += os.path.basename(im).replace(".fits", "") + ","
    imlist = imlist[:-1] + "]"

    wgts = sorted(
        glob(
            f"{datadir}/*.weights.{stoke.lower()}.*.beam*[00-35.fits]"
        )
    )
    weightlist = "["
    for wgt in wgts:
        weightlist += os.path.basename(wgt).replace(".fits", "") + ","
    weightlist = weightlist[:-1] + "]"

    parset_dir = datadir

    parset_file = f"{parset_dir}/linmos_{stoke}.in"
    parset = f"""linmos.names            = {imlist}
linmos.weights          = {weightlist}
linmos.imagetype        = fits
linmos.outname          = {ims[0][:ims[0].find('beam')]}linmos
linmos.outweight        = {wgts[0][:wgts[0].find('beam')]}linmos
linmos.weighttype       = Combined
linmos.weightstate      = Inherent
# Reference image for offsets
linmos.feeds.centre     = [{septab['FOOTPRINT_RA'][0]}, {septab['FOOTPRINT_DEC'][0]}]
linmos.feeds.spacing    = 1deg
# Beam offsets
"""
    for im in ims:
        basename = os.path.basename(im).replace('.fits', '')
        idx = basename.find('beam')
        beamno = int(basename[len('beam')+idx:len('beam')+idx+2])
        idx = septab['BEAM'] == beamno
        offset = f"linmos.feeds.{basename} = [{septab[idx]['DELTA_RA']},{septab[idx]['DELTA_DEC']}]\n"
        parset += offset
    with open(parset_file, "w") as f:
        f.write(parset)

    return parset_file


def genslurm(parsets, fieldname, cutdir, files, stokeslist, dryrun=True):
    """Generate SBATCH file

    Args:
        parsets (list): List containing path to all parsets
        fieldname (str): Name of RACS field
        cutdir (str): Directory containing cutouts
        files (list): List of directories inside cutouts
        dryrun (bool, optional): Generate SBATCH but do not submit. Defaults to True.
    """
    if dryrun:
        print("Doing a dryrun - no jobs will be submitted")
    slurm_dir = f"{cutdir}/slurmFiles"
    try:
        os.mkdir(slurm_dir)
    except FileExistsError:
        pass
    logdir = f"{slurm_dir}/outputs"
    try:
        os.mkdir(logdir)
    except FileExistsError:
        pass

    slurmbats = []
    for stoke in stokeslist:
        slurmbat = f"{slurm_dir}/{fieldname}_linmos_{stoke}.sbatch"
        slurm = f"""#!/bin/bash -l
#SBATCH --partition=workq
#SBATCH --clusters=galaxy
#SBATCH --exclude=nid00[010,200-202]
#SBATCH --account=askap
# No further constraints applied
# No reservation requested
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --time=24:00:00
#SBATCH --ntasks=20
#SBATCH --ntasks-per-node=20
#SBATCH --job-name={stoke}_{fieldname}_linmos
#SBATCH --export=NONE
#SBATCH --output={logdir}/slurm-%x.%j.out
#SBATCH --error={logdir}/slurm-%x.%j.err
log={logdir}/{stoke}_{fieldname}_linmos_$SLURM_JOB_ID.log

# Need to load the slurm module directly
module load slurm
# Ensure the default python module is loaded before askapsoft
module unload python
module load python
# Using user-defined askapsoft module
module use /group/askap/modulefiles
module load askapdata
module unload askapsoft
module load numpy
module load matplotlib
module load astropy
module load askapsoft
## MW experimental LINMOS
# module unload askapsoft
# module use /group/askap/wie017/modulefiles
# module load askapsoft/dev-omp
# Fixed linmos
#module load askapsoft/66f1e70
# Exit if we could not load askapsoft
#if [ "$ASKAPSOFT_RELEASE" == "" ]; then
#    echo "ERROR: \$ASKAPSOFT_RELEASE not available - could not load askapsoft module."
#    exit 1
#fi

NCORES=20
NPPN=20
mongo_ip=$(cat {logdir}/mongo_ip.txt)

task(){{
    dir={cutdir}/$1
    cd $dir
    linmos -c ${{dir}}/linmos_{stoke}.in >> "$log"
    #srun --export=ALL --ntasks=${{NCORES}} --ntasks-per-node=${{NPPN}} linmos-mpi -c ${{dir}}/linmos_{stoke}.in > "$log"
    ls ${{dir}}/*.cutout.sm.image.restored.{stoke.lower()}*.linmos.fits | xargs -I // mongo --host $mongo_ip --eval 'db.beams.findOneAndUpdate({{"Source_ID" : "'$1'"}}, {{"$set" :{{"beams.{fieldname}.{stoke.lower()}_file" : "'//'"}}}});' spiceracs >> "$log"
    echo $1
}}
"""
        # Get island dirs
        islands = []
        for file in files:
            base = os.path.basename(file)
            if base == "slurmFiles":
                continue
            else:
                islands.append(base)

        # Put them in a bash list
        island_list = ""
        for isl in islands:
            island_list += f"{isl} "

        cmd = f"""
islandList="{island_list}"
for island in $islandList; do
    ((i=i%NPPN)); ((i++==0)) && wait
    task $island &
done | tqdm --total {len(islands)} >> /dev/null
"""
        slurm += cmd
        with open(slurmbat, "w") as f:
            f.write(slurm)
        print(f"Written slurm file to {slurmbat}")
        slurmbats += [slurmbat]
    return slurmbats, slurm_dir, logdir


def yes_or_no(question):
    """Ask a yes or no question

    Args:
        question (str): The question

    Returns:
        bool: True for yes, False for no
    """
    while "Please answer 'y' or 'n'":
        reply = str(input(question+' (y/n): ')).lower().strip()
        if reply[:1] == 'y':
            return True
        if reply[:1] == 'n':
            return False


def main(args):
    """Main script
    """

    # Use ASKAPcli to get beam separations for PB correction
    field = args.field
    scriptdir = os.path.dirname(os.path.realpath(__file__))

    beamseps = gen_seps(field)

    stokeslist = args.stokeslist
    if stokeslist is None:
        stokeslist = ["I", "Q", "U", "V"]

    dryrun = args.dryrun

    if not dryrun:
        print('In the words of CA... check yoself before you wreck yoself!')
        dryrun = not yes_or_no(
            "Are you sure you want to submit jobs to the queue?")

    cutdir = args.cutdir
    if cutdir is not None:
        if cutdir[-1] == '/':
            cutdir = cutdir[:-1]

    files = sorted(glob(f"{cutdir}/*"))

    parfiles = []
    for file in tqdm(files, desc='Making parsets'):
        if os.path.basename(file) == 'slurmFiles':
            continue
        for stoke in stokeslist:
            parfile = genparset(field, stoke.capitalize(),
                                file, beamseps, prefix=args.prefix)
            parfiles.append(parfile)

    slurmbats, slurm_dir, logdir = genslurm(parfiles,
                                            field,
                                            cutdir,
                                            files,
                                            stokeslist)
    mongobat = gen_mongo(field, slurm_dir, logdir)
    gen_pipline(slurm_dir, field, mongobat, slurmbats, dryrun=dryrun)

    print('Done!')


def cli():
    """Command-line interface
    """
    import argparse

    # Help string to be shown using the -h option
    descStr = f"""
    Mosaic RACS beam cubes with linmos.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "field",
        metavar="field",
        type=str,
        help="RACS field to mosaic - e.g. 2132-50A."
    )

    parser.add_argument(
        "cutdir",
        metavar="cutdir",
        type=str,
        help="Directory containing the cutouts.",
    )

    parser.add_argument(
        "-d",
        "--dryrun",
        dest="dryrun",
        action="store_true",
        help="DON'T submit jobs (just make parsets) [False].",
    )

    parser.add_argument(
        "--prefix",
        metavar="prefix",
        type=str,
        default="",
        help="Prepend prefix to file.",
    )

    parser.add_argument(
        "-s",
        "--stokes",
        dest="stokeslist",
        nargs='+',
        type=str,
        help="List of Stokes parameters to image [ALL]",
    )

    args = parser.parse_args()

    main(args)


if __name__ == "__main__":
    cli()
