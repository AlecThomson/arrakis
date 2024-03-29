#!/bin/bash -l
#SBATCH --job-name=SPICE-CASDA
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/arrakis/spica/slurmLogs/casda_prep_slurm-%j.log
#SBATCH -o /group/askap/athomson/projects/arrakis/spica/slurmLogs/casda_prep_slurm-%j.log
#SBATCH --cluster=galaxy
#SBATCH --account=askap
#SBATCH --ntasks=1000
#SBATCH --ntasks-per-node=10
##SBATCH --time=0-00:45:00 # For cut
##SBATCH --time=0-00:10:00 # For test
#SBATCH --time=0-01:45:00 # For full

prep_type=full
conda activate spice
data_dir=/group/askap/athomson/projects/arrakis/DR1/full_spica
polcat=/group/askap/athomson/projects/arrakis/DR1/Arrakis.dr1.corrected.xml


cd $data_dir
srun -n $SLURM_NTASKS casda_prepare.py $data_dir $polcat $prep_type --convert-spectra --convert-cubes --convert-plots -v --mpi --batch_size 10_000
