#!/bin/bash -l

#SBATCH --job-name=racs_cuts
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --account=askap
#SBATCH -M galaxy
#SBATCH -p workq
#SBATCH --nodes=100
#SBATCH --ntasks=2000
#SBATCH --time=24:00:00

export OMP_NUM_THREADS=1

module load python/3.6.3
module load mpi4py
module unload python/3.6.3
conda activate py36

cd /group/askap/athomson/repos/spiceracs

numactl --interleave=all mongod --dbpath=database --bind_ip $(hostname -i) &

srun -n 2000 python spiceracs/cutout_rolling.py /group/askap/athomson/projects/RACS/VAST_2132-50A 2132-50A $(hostname -i) -v -vw --mpi

echo 'done'