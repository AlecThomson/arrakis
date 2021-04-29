#!/bin/bash -l

#SBATCH --job-name=processSPICE
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00

export OMP_NUM_THREADS=1

module load askapsoft
unset PYTHONPATH
source /home/$(whoami)/.bashrc
conda activate base


cd /group/askap/athomson/projects/spiceracs/testing_april_2021

# Assuing defaut port for Dask
# Forward the ports to the login nodes
port=9999
ssh -N -f -R $port:localhost:$port galaxy-1
ssh -N -f -R $port:localhost:$port galaxy-2

# Start MongoDB
host=$(hostname -i)
echo $host
database=/group/askap/athomson/projects/spiceracs/spiceracsDB
numactl --interleave=all mongod --dbpath=$database --bind_ip $host >> /dev/null &

processSPICE 1600-50A /home/athomson/athomson/projects/spiceracs/testing_april_2021/1600-50A $host --config test_config.txt