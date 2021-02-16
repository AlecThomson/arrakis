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

cd /group/askap/athomson/repos/spiceracs

# Assuing defaut port for Dask
# Forward the ports to the login nodes
port=8787
/usr/bin/ssh -N -f -R $port:localhost:$port galaxy-1
/usr/bin/ssh -N -f -R $port:localhost:$port galaxy-2

# Start MongoDB
host=$(hostname -i)
echo $host
numactl --interleave=all mongod --dbpath=database --bind_ip $host >> /dev/null &

python spiceracs/processSPICE.py 0918+06A /group/askap/athomson/projects/RACS/newpipelinetest/RACS_test4_1.05_0918+06A database/ $host -v -vw -p 0 -t -s 100 -o 5 -c -5 -m