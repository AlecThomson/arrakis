#!/bin/bash -l

#SBATCH --job-name=racs_convolve
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --account=askap
#SBATCH -M galaxy
#SBATCH -p workq
#SBATCH --ntasks=288
#SBATCH --time=24:00:00

export OMP_NUM_THREADS=1

module load slurm

module load convolve

cd /group/askap/athomson/projects/RACS/pipelinetest/RACS_test4_1.05_0918+06A


NCORES=288
NPPN=3

srun --pty --export=ALL -n 288 beamcon_3D.py image.restored*contcube*beam{00..35}.fits --mode total -v --bmaj 25 --bmin 25 --bpa 0 -c 25 --uselogs

echo 'done'