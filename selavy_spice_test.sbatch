#!/bin/bash -l
#SBATCH --partition=workq
#SBATCH --clusters=galaxy
#SBATCH --account=askap
#SBATCH --export=NONE
#SBATCH --ntasks=19
#SBATCH --ntasks-per-node=20
#SBATCH --time=08:00:00
#SBATCH --job-name=spice_test

module load slurm
module unload python
module load python
module use /group/askap/modulefiles
module load askapdata
module unload askapsoft
module load askapsoft/0.24.6
module unload askappipeline
module load askappipeline/0.24.6

srun --export=ALL --ntasks=19 --ntasks-per-node=20 selavy -c selavy_spice_test.in >> selavy_${{SLURM_JOB_ID}}.log
