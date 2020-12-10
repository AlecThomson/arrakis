#!/bin/bash -l

#SBATCH --job-name=SPICE-RMsynth
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=1
#SBATCH --time=24:00:00

cd /group/askap/athomson/repos/spiceracs

module load python
module load mpi4py
module unload python
conda activate py36

numactl --interleave=all mongod --dbpath=database --bind_ip $(hostname -i) >> /dev/null &

srun -n 1 -c 20 python spiceracs/rmsynth_oncuts.py 2132-50A ~/athomson/projects/RACS/VAST_2132-50A/ $(hostname -i) -v --ncores 20 -m -t -l 1000 -s 100

echo 'Done'