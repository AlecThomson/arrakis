#!/bin/bash -l

#SBATCH --job-name=SPICE-CASDA
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.err
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.out
#SBATCH --ntasks=500
#SBATCH --ntasks-per-node=10

#SBATCH --time=1-00:00:00
#SBATCH --cluster=galaxy
#SBATCH --account=askap

##SBATCH --cluster=magnus
##SBATCH --account=ja3
##SBATCH -p debugq

conda activate spice

data_dir=/group/ja3/athomson/full_spica
polcat=/group/ja3/athomson/spice-racs.dr1.corrected.xml

cd $data_dir
srun casda_prepare.py $data_dir $polcat --update-cubes --convert-spectra --convert-plots -v --mpi --batch_size $SLURM_NTASKS