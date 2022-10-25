#!/bin/bash -l

#SBATCH --job-name=SPICE-CASDA
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/casda_pre_prep_slurm-%j.log
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/casda_pre_prep_slurm-%j.log

#SBATCH --cluster=zeus
#SBATCH --account=askap
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=2-00:00:00
#SBATCH --partition=copyq

data_dir=/askapbuffer/processing/tho822/spice-racs/DR1/full_spica/cutouts
polcat=/askapbuffer/processing/tho822/spice-racs/DR1/spice-racs.dr1.corrected.xml

data_dir_s=/group/askap/athomson/projects/spiceracs/DR1/full_spica/cutouts
polcat_s=/group/askap/athomson/projects/spiceracs/DR1/spice-racs.dr1.corrected.xml


module load rclone

rclone sync -P --transfers 10 --checkers 10 $data_dir $data_dir_s
rclone sync -P --transfers 10 --checkers 10 $polcat $polcat_s