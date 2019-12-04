#!/bin/bash -l
#SBATCH --partition=workq
#SBATCH --clusters=galaxy
#SBATCH --account=askap
#SBATCH --export=NONE
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=spicemoment
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alec.thomson@anu.edu.au
module load slurm

cd /group/askap/athomson/

spicemoments /group/askap/mcc381/RACS/8585/alec/ /group/askap/mcc381/RACS/8585/RACS_test4_1.05_1049-31A/selavy-cont-image.i.SB8583.cont.RACS_test4_1.05_1049-31A.linmos.taylor.0.restored/ /group/askap/mcc381/RACS/8585/RACS_test4_1.05_1049-31A/ /group/askap/athomson/ -v --zero 1>prog_moment.out 2>&1
