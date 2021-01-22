#!/bin/bash -l

#SBATCH --job-name=SPICE-RMsynth
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=5
#SBATCH --time=24:00:00

cd /group/askap/athomson/repos/spiceracs

module load python
module load mpi4py
module unload python
conda activate py36

host=$(hostname -i)

echo 'Mongo running on '$host

numactl --interleave=all mongod --dbpath=database --bind_ip $host >> /dev/null &

srun -n 100 python spiceracs/rmsynth_oncuts.py 0918+06A /home/athomson/athomson/projects/RACS/RACS_test4_1.05_0918+06A $host -v -m -t -l 1000 -s 100 --mpi

echo 'Done'