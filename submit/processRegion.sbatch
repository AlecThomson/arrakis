#!/bin/bash -l

#SBATCH --job-name=processSPICE
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.err
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/slurm-%j.out
#SBATCH --ntasks=480

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

srun -n 480 --export=ALL python spiceracs/process_region.py --merge_name test --fields 1212-12A 1237-12A --datadirs /group/ja3/athomson/spica/12509/RACS_test4_1.05_1212-12A/ /group/ja3/athomson/spica/12509/RACS_test4_1.05_1237-12A/ --output_dir /group/ja3/athomson --host 146.118.68.63 --username admin --password 9skBMPshVXk4JWi2HC -v --use_mpi --tt0 /group/askap/athomson/projects/RACS/CI0_mosaic_1.0/spica.tt0.fits --tt1 /group/askap/athomson/projects/RACS/CI1_mosaic_1.0/spica.tt1.fits --port_forward magnus-1 magnus-2