#!/bin/bash -l

#SBATCH --job-name=processSPICE
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=7:00:00
#SBATCH --cluster=magnus
#SBATCH --account=ja3
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.err
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.out

export OMP_NUM_THREADS=1

# Module sorting
module unload askapsoft
module load askapsoft
unset PYTHONPATH
source /home/$(whoami)/.bashrc
conda activate spice


# Image dirctory
cd /scratch/ja3/athomson/spica

# Start MongoDB
host=$(hostname -i)
echo $host
database=/group/askap/athomson/projects/spiceracs/spica/spicaDB
numactl --interleave=all mongod --dbpath=$database --bind_ip $host >> /dev/null &

field=1213-25A
cal_sbid=8577

processSPICE $field /scratch/ja3/athomson/spica/$cal_sbid/RACS_test4_1.05_$field $host --config /group/askap/athomson/projects/spiceracs/spica/spica_config.txt --savePlots