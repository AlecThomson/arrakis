#!/bin/bash -l

#SBATCH --job-name=racs_smooth
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --account=askap
#SBATCH -M galaxy
#SBATCH -p workq
#SBATCH --nodes=15
#SBATCH --ntasks=300
#SBATCH --time=24:00:00

export OMP_NUM_THREADS=1

module load python/3.6.3
module load numpy
module load astropy
module load scipy
module load mpi4py
conda activate py36
module load pandas

cd /group/askap/athomson/projects/RACS/racs_smooth

numactl --interleave=all mongod --dbpath database/ &

srun -n 288 python /group/askap/athomson/repos/spiceracs/spiceracs/cutout.py ../RACS_test4_1.05_1145-31A ../RACS_test4_1.05_1145-31A/selavy-cont-image.i.SB8593.cont.RACS_test4_1.05_1145-31A.linmos.taylor.0.restored/ ../RACS_test4_1.05_1145-31A/selavy-cont-image.i.SB8593.cont.RACS_test4_1.05_1145-31A.linmos.taylor.0.restored/snrMap.image.i.SB8593.cont.RACS_test4_1.05_1145-31A.linmos.taylor.0.restored.fits . 5 -v --mpi