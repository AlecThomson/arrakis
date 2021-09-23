#!/bin/bash -l

#SBATCH --job-name=processSPICE
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=24
#SBATCH --nodes=20
#SBATCH --ntasks=480
#SBATCH --time=01:00:00
#SBATCH --cluster=magnus
#SBATCH --account=ja3
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.err
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.out

export OMP_NUM_THREADS=1

# Module sorting
# module unload askapsoft
# module load askapsoft
# unset PYTHONPATH
# source /home/$(whoami)/.bashrc
conda activate spice
module load singularity
export SINGULARITY_BINDPATH=$(pwd)


# Image dirctory
cd /scratch/ja3/athomson/spica

field=1416+00A
cal_sbid=13615
tt0_dir=/group/askap/athomson/projects/RACS/CI0_mosaic_1.0
tt1_dir=/group/askap/athomson/projects/RACS/CI1_mosaic_1.0

srun -n 480 --export=ALL processSPICE $field /scratch/ja3/athomson/spica/$cal_sbid/RACS_test4_1.05_$field --config /group/askap/athomson/projects/spiceracs/spica/spica_config.txt --savePlots --tt0 $tt0_dir/RACS_test4_1.05_$field.fits --tt1 $tt1_dir/RACS_test4_1.05_$field.fits --use_mpi --skip_cutout