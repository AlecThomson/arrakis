#!/bin/bash -l

#SBATCH --job-name=processSPICE
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/arrakis/spica/slurmLogs/slurm-%j.err
#SBATCH -o /group/askap/athomson/projects/arrakis/spica/slurmLogs/slurm-%j.out
#SBATCH --ntasks=480

#SBATCH --ntasks-per-node=20
#SBATCH --time=1-00:00:00
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
export SINGULARITY_BINDPATH=$(pwd),/group,/askapbuffer

field=1335-25A
sbid=`find_sbid.py $field --science`
tt0_dir=/group/askap/athomson/projects/RACS/CI0_mosaic_1.0
tt1_dir=/group/askap/athomson/projects/RACS/CI1_mosaic_1.0
data_dir=/askapbuffer/processing/len067/arrakis
config=/group/askap/athomson/projects/arrakis/spica/spica_config.txt

# Image dirctory
cd $data_dir

# Make a copy of this sbatch file for posterity
sedstr="s/sbatch/${field}.${SLURM_JOB_ID}\.sbatch/g"
slurmdir="/group/askap/athomson/projects/arrakis/spica/slurmFiles"
currentdir="/group/askap/athomson/repos/arrakis/submit"
sedstr2="s+${currentdir}+${slurmdir}+g"
thisfile="/group/askap/athomson/repos/arrakis/submit/process_reSPICE.sh"
cp "$thisfile" "$(echo "$thisfile" | sed -e "$sedstr" | sed -e "$sedstr2")"

srun -n 480 --export=ALL processSPICE $field $data_dir/$sbid/RACS_test4_1.05_$field --config $config --savePlots --tt0 $tt0_dir/RACS_test4_1.05_$field.fits --tt1 $tt1_dir/RACS_test4_1.05_$field.fits --use_mpi --port_forward galaxy-1 galaxy-2 -v