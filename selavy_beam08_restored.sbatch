#!/bin/bash -l
#SBATCH --partition=workq
#SBATCH --clusters=magnus
#SBATCH --account=ja3
#SBATCH --mail-user=andrew.obrien@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --time=05:00:00
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=20
#SBATCH --job-name=selavyContR_B08
#SBATCH --export=NONE

module load slurm
module unload python
module load python
module use /group/askap/modulefiles
module load askapdata
module unload askapsoft
module load askapsoft/0.24.4
module unload askappipeline
module load askappipeline/0.24.4

NCORES=2
NPPN=20
parset="selavy_beam08_restored.in"
log="selavy_beam08_restored_${SLURM_JOB_ID}.log"
srun --export=ALL --ntasks=${NCORES} --ntasks-per-node=${NPPN} selavy -c "$parset" >> "$log"
