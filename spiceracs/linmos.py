#!/usr/bin/env python
import os
import numpy as np
from astropy.table import Table
from glob import glob
from astropy.coordinates import SkyCoord
import astropy.units as u
import subprocess
import shlex


def genparset(field, stoke, datadir, prefix=""):

    ims = sorted(glob(f"{datadir}/{prefix}image.restored.{stoke.lower()}.*.contcube.VAST_{field}.beam*.fits"))
    if len(ims)==0:
        print(f"{datadir}/{prefix}image.restored.{stoke.lower()}.*.contcube.VAST_{field}.beam*.fits")
        raise Exception('No files found. Have you run imaging? Check your prefix?')
    imlist = "["
    for im in ims:
        imlist += os.path.basename(im).replace(".fits", "") + ","
    imlist = imlist[:-1] + "]"

    wgts = sorted(glob(
        f"{datadir}/weights.{stoke.lower()}.*.contcube.VAST_{field}.beam*.fits"
    ))
    weightlist = "["
    for wgt in wgts:
        weightlist += os.path.basename(wgt).replace(".fits", "") + ","
    weightlist = weightlist[:-1] + "]"

    parset_dir = f"{scriptdir}/parsets/VAST_{field}"
    try:
        os.mkdir(parset_dir)
    except FileExistsError:
        pass

    parset_file = f"{parset_dir}/science_contcube_linmos_F04_{stoke}.in"
    parset = f"""linmos.names            = {imlist}
linmos.weights          = {weightlist}
linmos.imagetype        = fits
linmos.outname          = {ims[0][:ims[0].find('beam')]}linmos
linmos.outweight        = {wgts[0][:wgts[0].find('beam')]}linmos
linmos.weighttype       = Combined
linmos.weightstate      = Inherent
# Reference image for offsets
linmos.feeds.centre     = [10h17m49.958, -27.49.24.30]
linmos.feeds.spacing    = 1deg
# Beam offsets
linmos.feeds.image.i.SB10269.cont.Hydra_1B.beam13.taylor.0.restored = [-1.344,-1.356]
linmos.feeds.image.i.SB10269.cont.Hydra_1B.beam32.taylor.0.restored = [-2.244,-1.359]
    """
    with open(parset_file, "w") as f:
        f.write(parset)

    return parset_file

def genslurm(dryrun=True):
    if dryrun:
        print("Doing a dryrun - no jobs will be submitted")
    slurm_dir = f"{scriptdir}/slurmFiles/VAST_{field}"
    try:
        os.mkdir(slurm_dir)
    except FileExistsError:
        pass
    logdir = f"{scriptdir}/outputs"
    try:
        os.mkdir(logdir)
    except FileExistsError:
        pass
    slurmbat = f"{slurm_dir}/science_contcube_linmos_F04_{stoke}.sbatch"
    slurm = f"""#!/bin/bash -l
#SBATCH --partition=workq
#SBATCH --clusters=galaxy
#SBATCH --exclude=nid00[010,200-202]
#SBATCH --account=askap
# No further constraints applied
# No reservation requested
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --time=10:00:00
#SBATCH --ntasks=288
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=linmosCCrestored_F04_{stoke}
#SBATCH --export=NONE
#SBATCH --output={logdir}/slurm-%x.%j.out
#SBATCH --error={logdir}/slurm-%x.%j.err
log={logdir}/RACS_{field}_science_linmosCC_F04_{stoke}_$SLURM_JOB_ID.log

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

cd {datadir}

NCORES=288
NPPN=4
srun --export=ALL --ntasks=$NCORES --ntasks-per-node=$NPPN linmos-mpi -c {parset_file} > "$log"
err=$?
    """
    with open(slurmbat, "w") as f:
        f.write(slurm)

    if not dryrun:
        print(f"Submitting {slurmbat}")
        command = f"sbatch {slurmbat}"
        command = shlex.split(command)
        proc = subprocess.run(
            command, capture_output=True, encoding="utf-8", check=True
        )


def yes_or_no(question):
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
    SBID = args.SBID
    field = args.field

    # Get the SBID info
    stdout = subprocess.getoutput(f'schedblock info -p {SBID}') 
    info = stdout.splitlines()
    for line in info[7:]: 
        temp = line.split() 
        if temp[-1] == '=': 
            info_dict.update( 
                { 
                    temp[0]: None 
                } 
            ) 
        else: 
            info_dict.update( 
                { 
                    temp[0]: temp[-1]
                } 
            ) 

    keys, vals = [], []
    for line in info[7:]: 
        line_lst = line.split() 
        val = line_lst[-1] 
        key = line_lst[0].split('.')
        keys.append(key) 
        vals.append(val)

    srcs = [] 
    for key, val in zip(keys, vals):
        if key[0] == 'craft': 
            continue 
        elif key[1] == 'target': 
            srcs.append(key[2])

    # Store data in dict stucture
    srcs = sorted(np.unique(srcs))
    info_dict = {}
    for src in srcs: 
        info_dict.update({src: {} })

    info_dict['src%d']['footprint'] = {}

    for key, val in zip(keys, vals): 
        if key[0] == 'craft':  
            continue  
        elif key[1] == 'target':  
            src = key[2]
            if len(key[3:]) == 1:
                    info_dict[src].update(
                        {
                            key[3]: val
                        }
                    )
            elif key[3] == 'footprint':
                print(key[3], key[4], val)
                info_dict[src]['footprint'].update(
                    {
                        key[4]: val
                    }
                )

    # Get the key direction and rotation info
    for src in info_dict.keys(): 
        if info_dict[src]['field_name'][:-len(field)] == field: 
            pol_lst = info_dict[src]['pol_axis'].split(',') 
            rot_str = pol_lst[1] 
            idx = rot_str.find(']') 
            pol_ax_rot = float(rot_str[0:idx]) 
            foot_rot = float(info_dict['src%d']['footprint']['rotation'])
            rotation = pol_ax_rot+foot_rot
            direction = info_dict[src]['field_direction'][len('['):-len(']')].split(',')                         

    # Get the footprint seperations for each beam
    cmd = f"footprint calculate -d {direction[0]},{direction[0]} -r {rotation} \
        {info_dict['src%d']['footprint']['name']} -p {info_dict['src%d']['footprint']['pitch']}"
    stdout = subprocess.getoutput(cmd)
    beam_sep = {}
    for beam, line in enumerate(stdout.splitlines()):
        start = line.find('(')
        stop = line.find(')')
        beam_sep.update(
            {
                beam: line[start+1:stop].split()
            }
        )


    stokeslist = args.stokeslist
    if stokeslist is None:
        stokeslist = ["I", "Q", "U", "V"]

    dryrun = args.dryrun

    if not dryrun:
        print('In the words of CA... check yoself before you wreck yoself!')
        dryrun = not yes_or_no("Are you sure you want to submit jobs to the queue?")

    cutdir = args.cutdir
    if cutdir is not None:
        if cutdir[-1] == '/':
            cutdir = cutdir[:-1]

    files = sorted(glob(f"{cutdir}/*"))    

    parfiles = []
    for file in files:
        for stoke in stokeslist:
            parfile = genparset(field, stoke.capitalize(), file, prefix=args.prefix)
            parfiles.append(parfile)


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
        "SBID", 
        metavar="SBID", 
        type=str, 
        help="RACS SBID corresponding to field."
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
