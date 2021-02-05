#!/bin/bash -l

#SBATCH --job-name=processSPICE
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00


cd /group/askap/athomson/repos/spiceracs

python spiceracs/processSPICE.py /group/askap/athomson/projects/RACS/pipelinetest/RACS_test4_1.05_0918+06A 0918+06A database/ -v -vw