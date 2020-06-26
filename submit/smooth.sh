#!/bin/bash -l

#SBATCH --job-name=racs_smooth
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --account=askap
#SBATCH -M galaxy
#SBATCH -p workq
#SBATCH --nodes=50
#SBATCH --ntasks=1000
#SBATCH --time=24:00:00

export OMP_NUM_THREADS=1

module load python/3.6.3
module load mpi4py
module unload python/3.6.3
conda activate py36

cd /group/askap/athomson/repos/spiceracs

srun -n 1000 python spiceracs/beamcon.py /group/askap/athomson/projects/RACS/VAST_2132-50A/cutouts/ /group/askap/athomson/projects/RACS/VAST_2132-50A/ -v --bmaj 25 --bmin 25 --bpa 0 -c 25 -vw --mpi

echo 'done'