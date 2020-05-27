#!/bin/bash -l

#SBATCH --job-name=racs_cuts
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --account=askaprt
#SBATCH --partition=highmemq
#SBATCH --nodelist=z140
#SBATCH --ntasks=16
#SBATCH --time=24:00:00

export OMP_NUM_THREADS=1

source /group/askap/athomson/miniconda3/bin/activate
conda activate py36
module load python/3.6.3
module load numpy
module load astropy
module load scipy
module unload sandybridge
module load mpi4py
module load pandas
conda activate py36

TEMPDIR=/dev/shm
DATADIR=/group/askap/athomson/projects/RACS/VAST_2321-50A
POLLST="i q u v"

cd $DATADIR

numactl --interleave=all /group/askap/athomson/miniconda3/bin/mongod --dbpath database/ &

echo Copying cubes to $TEMPDIR
for pol in $POLLST; do
    time cp sm.image.restored.$pol.SB11628.contcube.VAST_2132-50A.linmos.fits $TEMPDIR/image.restored.$pol.SB11628.contcube.VAST_2132-50A.linmos.fits
done

echo Starting cutout
srun -n 1 -c 14 python /group/askap/athomson/repos/spiceracs/spiceracs/cutout.py $TEMPDIR $DATADIR/selavy/ $DATADIR/selavy/image.i.SB8644.cont.RACS_test4_1.05_2132-50A.linmos.taylor.0.restored.sm.fits $TEMPDIR 3 -v --ncores 14 -m


echo Copying cutouts to $DATADIR
time cp -r $TEMPDIR/cutouts $DATADIR


/group/askap/athomson/miniconda3/bin/mongod --shutdown --dbpath database/

echo 'Done'