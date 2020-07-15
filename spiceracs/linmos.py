#!/usr/bin/env python

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
            f"{datadir}/sm.*.cutout.image.restored.{stoke.lower()}.*.beam*.fits"
        )
    )
    if len(ims) == 0:
        print(
            f"{datadir}/sm.*.cutout.image.restored.{stoke.lower()}.*.beam*.fits"
        )
        raise Exception(
            'No files found. Have you run imaging? Check your prefix?')
    imlist = "["
    for im in ims:
        imlist += os.path.basename(im).replace(".fits", "") + ","
    imlist = imlist[:-1] + "]"

    wgts = sorted(
        glob(
            f"{datadir}/*.weights.{stoke.lower()}.*.beam*.fits"
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
        beamno = int(basename[len('beam')+idx:])
        idx = septab['BEAM'] == beamno
        offset = f"linmos.feeds.{basename} = [{septab[idx]['DELTA_RA'][0]},{septab[idx]['DELTA_DEC'][0]}]\n"
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
    slurmbat = f"{slurm_dir}/{fieldname}_linmos.sbatch"
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
#SBATCH --ntasks=288
#SBATCH --ntasks-per-node=20
#SBATCH --job-name={fieldname}_linmos
#SBATCH --export=NONE
#SBATCH --output={logdir}/slurm-%x.%j.out
#SBATCH --error={logdir}/slurm-%x.%j.err
log={logdir}/{fieldname}_linmos_$SLURM_JOB_ID.log

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
# Exit if we could not load askapsoft
if [ "$ASKAPSOFT_RELEASE" == "" ]; then
    echo "ERROR: \$ASKAPSOFT_RELEASE not available - could not load askapsoft module."
    exit 1
fi

NCORES=288
NPPN=20

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
    dir={cutdir}/${{island}}
    cd $dir
    srun --export=ALL --ntasks=$NCORES --ntasks-per-node=$NPPN linmos-mpi -c ${{dir}}/linmos_I.in >> "$log"
    srun --export=ALL --ntasks=$NCORES --ntasks-per-node=$NPPN linmos-mpi -c ${{dir}}/linmos_Q.in >> "$log"
    srun --export=ALL --ntasks=$NCORES --ntasks-per-node=$NPPN linmos-mpi -c ${{dir}}/linmos_U.in >> "$log"
    srun --export=ALL --ntasks=$NCORES --ntasks-per-node=$NPPN linmos-mpi -c ${{dir}}/linmos_V.in >> "$log"
done | tqdm --total {len(islands)} >> /dev/null
"""
    slurm += cmd
    with open(slurmbat, "w") as f:
        f.write(slurm)
    print(f"Written slurm file to {slurmbat}")
    if not dryrun:
        print(f"Submitting {slurmbat}")
        command = f"sbatch {slurmbat}"
        command = shlex.split(command)
        proc = subprocess.run(
            command, capture_output=True, encoding="utf-8", check=True
        )


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
    sepfile = f"{scriptdir}/../askap_surveys/RACS/admin/epoch_0/beam_sep-RACS_test4_1.05_{field}.csv"

    beamseps = Table.read(sepfile)

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
    
    genslurm(parfiles, field, cutdir, files, stokeslist, dryrun=dryrun)

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
