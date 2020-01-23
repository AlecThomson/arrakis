#!/bin/bash
#SBATCH --job-name=SPICE-CUT
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=1
#SBATCH --mem=256g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alec.thomson@csiro.au

# Application specific commands:
cd /flush3/tho822/cubes/8594

mpirun -np 20 spicecutout mosaics/ selavy/selavy-cont-image.i.SB8593.cont.RACS_test4_1.05_1145-31A.linmos.taylor.0.restored/ selavy/selavy-cont-image.i.SB8593.cont.RACS_test4_1.05_1145-31A.linmos.taylor.0.restored/componentMap_image.i.SB8593.cont.RACS_test4_1.05_1145-31A.linmos.taylor.0.restored.fits ./ 2 -v -m --mpi 1>prog_cut.out 2>&1
