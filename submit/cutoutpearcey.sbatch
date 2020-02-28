#!/bin/bash -l

#SBATCH --job-name=cutout
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=15
#SBATCH --time=24:00:00


module load python/3.6.1
module load mongodb
module load openmpi

cd /flush3/tho822/common_racs/

numactl --interleave=all mongod --dbpath database/ &

mpirun -n 288 /apps/python/3.6.1/bin/python /home/tho822/repos/spiceracs/spiceracs/cutout.py mosaics selavy-cont-image.i.SB8593.cont.RACS_test4_1.05_1145-31A.linmos.taylor.0.restored/ selavy-cont-image.i.SB8593.cont.RACS_test4_1.05_1145-31A.linmos.taylor.0.restored/snrMap.image.i.SB8593.cont.RACS_test4_1.05_1145-31A.linmos.taylor.0.restored.fits . 5 -v --mpi