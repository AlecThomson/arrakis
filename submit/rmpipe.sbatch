#!/bin/bash -l

#SBATCH --job-name=SPICE-RM
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=18
#SBATCH --nodes=50
#SBATCH --time=24:00:00

cd /group/askap/athomson/repos/spiceracs

module load python
module load mpi4py
module unload python
conda activate py36

host=$(hostname -i)

echo 'Mongo running on '$host

numactl --interleave=all mongod --dbpath=database --bind_ip $host >> /dev/null &

echo 'Starting synth'

srun -n 900 python spiceracs/rmsynth_oncuts.py 2132-50A ~/athomson/projects/RACS/VAST_2132-50A/ $host -v -m -t -l 900 -s 100 --mpi -rmv

echo 'Starting CLEAN'

srun -n 900 python spiceracs/rmclean_oncuts.py 2132-50A ~/athomson/projects/RACS/VAST_2132-50A/ $host -v -m --mpi -rmv

echo 'Done'