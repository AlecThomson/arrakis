#!/bin/bash -l

#SBATCH --job-name=SPICE-MERGE
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/process_spice-%j.log
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/process_spice-%j.log
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
cd /askapbuffer/processing/tho822/spice-racs/DR1
conda activate spice
module load singularity
export SINGULARITY_BINDPATH=$(pwd),/group


srun --export=ALL spice_region --dask_config /group/askap/athomson/repos/spiceracs/configs/galaxy.yaml --config /group/askap/athomson/projects/spiceracs/spica/spica_full_region_config.txt --use_mpi --own_fit --outfile spice-racs.dr1.xml --skip_merge --output_dir /askapbuffer/processing/tho822/spice-racs/DR1

spice_script fix_dr1_cat spice-racs.dr1.xml