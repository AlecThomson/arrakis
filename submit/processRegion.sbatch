#!/bin/bash -l

#SBATCH --job-name=testSPICE
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.err
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.out
#SBATCH --ntasks=500

#SBATCH --ntasks-per-node=10
#SBATCH --time=12:00:00
#SBATCH --cluster=galaxy
#SBATCH --account=askap


export OMP_NUM_THREADS=1

# Module sorting
# module unload askapsoft
# module load askapsoft
# unset PYTHONPATH
# source /home/$(whoami)/.bashrc
cd /group/ja3/athomson
conda activate spice
module load singularity
export SINGULARITY_BINDPATH=$(pwd),/group


srun --export=ALL python /group/askap/athomson/repos/spiceracs/spiceracs/process_region.py --config /group/askap/athomson/projects/spiceracs/spica/spica_full_region_config.txt --use_mpi --skip_merge --skip_rmsynth --debugger