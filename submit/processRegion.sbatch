#!/bin/bash -l

#SBATCH --job-name=SPICE-MERGE
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/process_spice-%j.log
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/process_spice-%j.log
#SBATCH --ntasks=1000
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
cd /group/askap/athomson/projects/spiceracs/DR1
conda activate spice
module load singularity
export SINGULARITY_BINDPATH=$(pwd),/group
outfile=spice-racs.dr1.multiflag.test.xml

srun --export=ALL spice_region --dask_config /group/askap/athomson/repos/spiceracs/spiceracs/configs/galaxy.yaml --config /group/askap/athomson/projects/spiceracs/spica/spica_full_region_config.txt --use_mpi --own_fit --outfile $outfile --skip_merge --output_dir /group/askap/athomson/projects/spiceracs/DR1 --polyOrd -1

fix_dr1_cat.py $outfile