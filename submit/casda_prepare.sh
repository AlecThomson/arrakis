#!/bin/bash -l

#SBATCH --job-name=SPICE-CASDA
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.err
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.out

##SBATCH --ntasks=500
##SBATCH --ntasks-per-node=10
##SBATCH --time=0-03:00:00
##SBATCH --cluster=galaxy
##SBATCH --account=askap

##SBATCH --cluster=magnus
##SBATCH --account=ja3
##SBATCH --ntasks=48
##SBATCH --ntasks-per-node=12
##SBATCH --time=1-00:00:00

#SBATCH --cluster=zeus
#SBATCH --account=askap
#SBATCH --ntasks=48
#SBATCH --time=0-12:00:00
#SBATCH --partition=highmemq

# conda activate spice
module load intel-mpi
conda activate spice-zeus

data_dir=/group/ja3/athomson/full_spica
polcat=/group/ja3/athomson/spice-racs.dr1.corrected.xml

cd $data_dir
# srun -n $SLURM_NTASKS casda_prepare.py $data_dir $polcat --update-cubes --convert-spectra --convert-plots -v --mpi --batch_size 10_000
# srun -n $SLURM_NTASKS casda_prepare.py $data_dir $polcat --convert-plots -v --mpi --batch_size 10_000
# srun -n $SLURM_NTASKS casda_prepare.py $data_dir $polcat --convert-plots --update-cubes --convert-spectra -v --mpi --batch_size 10_000 --interface ib0 --outdir /scratch/ja3/athomson/spica
# srun -n $SLURM_NTASKS casda_prepare.py $data_dir $polcat --convert-plots --convert-spectra -v --mpi --batch_size 10_000 --interface ib0 --outdir /scratch/ja3/athomson/spica
srun -n $SLURM_NTASKS casda_prepare.py $data_dir $polcat --convert-spectra -v --mpi --batch_size 10_000 --interface ib0 --outdir /scratch/ja3/athomson/spica
# srun -n 20 casda_prepare.py $data_dir $polcat --convert-spectra --convert-plots -v --mpi --batch_size 1000
# srun -n 20 casda_prepare.py $data_dir $polcat --convert-spectra -v --mpi --batch_size 10_000