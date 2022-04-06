#!/bin/bash -l

#SBATCH --job-name=SPICE-MERGE
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.err
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.out
#SBATCH --ntasks=500
#SBATCH --requeue
#SBATCH --ntasks-per-node=10
#SBATCH --time=12:00:00
#SBATCH --cluster=galaxy
#SBATCH --account=askap
##SBATCH --dependency=afterok:13676803


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


srun --export=ALL spice_region --config /group/askap/athomson/projects/spiceracs/spica/spica_full_region_config.txt --use_mpi --skip_merge