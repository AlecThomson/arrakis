#!/bin/bash -l

#SBATCH --job-name=processSPICE
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.err
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.out
#SBATCH --ntasks=480

##SBATCH --ntasks-per-node=24
##SBATCH --time=06:00:00
##SBATCH --cluster=magnus
##SBATCH --account=ja3

#SBATCH --ntasks-per-node=20
#SBATCH --time=01:00:00
#SBATCH --cluster=galaxy
#SBATCH --account=askap


export OMP_NUM_THREADS=1

# Module sorting
# module unload askapsoft
# module load askapsoft
# unset PYTHONPATH
# source /home/$(whoami)/.bashrc
conda activate spice
module load singularity
export SINGULARITY_BINDPATH=$(pwd),/group

field=1213-25A
cal_sbid=`find_sbid.py $field --cal`
tt0_dir=/group/askap/athomson/projects/RACS/CI0_mosaic_1.0
tt1_dir=/group/askap/athomson/projects/RACS/CI1_mosaic_1.0
data_dir=/scratch/ja3/athomson/spica
# data_dir=/group/ja3/athomson/spica
config=/group/askap/athomson/projects/spiceracs/spica/spica_config.txt

# Image dirctory
cd $data_dir

# Make a copy of this sbatch file for posterity
sedstr="s/sbatch/${field}.${SLURM_JOB_ID}\.sbatch/g"
slurmdir="/group/askap/athomson/projects/spiceracs/spica/slurmFiles"
currentdir="/group/askap/athomson/repos/spiceracs/submit"
sedstr2="s+${currentdir}+${slurmdir}+g"
thisfile="/group/askap/athomson/repos/spiceracs/submit/processSPICE.sbatch"
cp "$thisfile" "$(echo "$thisfile" | sed -e "$sedstr" | sed -e "$sedstr2")"

srun -n 480 --export=ALL processSPICE $field $data_dir/$cal_sbid/RACS_test4_1.05_$field --config $config --savePlots --tt0 $tt0_dir/RACS_test4_1.05_$field.fits --tt1 $tt1_dir/RACS_test4_1.05_$field.fits --use_mpi