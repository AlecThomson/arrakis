#!/bin/bash -l

#SBATCH --job-name=SPICE-CASDA
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/casda_prep_slurm-%j.log
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/casda_prep_slurm-%j.log

#SBATCH --cluster=galaxy
#SBATCH --account=askap
#SBATCH --ntasks=100
#SBATCH --time=2-00:00:00

# conda activate spice
conda activate spice

data_dir=/group/askap/athomson/projects/spiceracs/DR1/full_spica
polcat=/group/askap/athomson/projects/spiceracs/DR1/spice-racs.dr1.corrected.xml

cd $data_dir
srun -n $SLURM_NTASKS casda_prepare.py $data_dir $polcat --convert-spectra --convert-cubes --convert-plots -v --mpi --batch_size 10_000