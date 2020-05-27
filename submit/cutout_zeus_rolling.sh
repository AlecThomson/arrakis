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
DATADIR=/group/askap/athomson/projects/RACS/VAST_2132-50A
POLLST="i q u"

cd /group/askap/athomson/repos/spiceracs

numactl --interleave=all mongod --dbpath database/ &

cd $DATADIR

echo Copying cubes to $TEMPDIR
for pol in $POLLST; do
    time cp image.restored.$pol.SB11628.contcube.VAST_2132-50A.beam*.fits $TEMPDIR/
    time cp weights.$pol.SB11628.contcube.VAST_2132-50A.beam*.fits $TEMPDIR/
done

echo Starting cutout
srun -n 14 python /group/askap/athomson/repos/spiceracs/spiceracs/cutout_rolling.py $TEMPDIR/ 2132-50A -v -vw --mpi

echo Copying cutouts to $DATADIR
time cp -r $TEMPDIR/cutouts $DATADIR

cd /group/askap/athomson/repos/spiceracs
/group/askap/athomson/miniconda3/bin/mongod --shutdown --dbpath database/

echo 'Done'