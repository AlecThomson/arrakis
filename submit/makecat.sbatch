#!/bin/bash -l

#SBATCH --job-name=SPICE-RMCLEAN
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00

cd /group/askap/athomson/repos/spiceracs

module load python
module load mpi4py
module unload python
conda activate py36

host=$(hostname -i)

echo 'Mongo running on '$host

numactl --interleave=all mongod --dbpath=database --bind_ip $host >> /dev/null &

python spiceracs/makecat.py 0918+06A . $host -w testcat.fits -f fits