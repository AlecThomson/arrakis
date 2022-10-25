#!/bin/bash -l

#SBATCH --job-name=SPICE-CASDA
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/casda_pre_prep_slurm-%j.log
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/casda_pre_prep_slurm-%j.log

#SBATCH --cluster=zeus
#SBATCH --account=askap
#SBATCH --ntasks=20
#SBATCH --time=1-00:00:00
#SBATCH --partition=copyq

data_dir=/askapbuffer/processing/tho822/spice-racs/DR1/full_spica/cutouts
polcat=/askapbuffer/processing/tho822/spice-racs/DR1/spice-racs.dr1.corrected.xml

data_dir_s=/group/askap/athomson/projects/spiceracs/DR1/full_spica/cutouts
polcat_s=/group/askap/athomson/projects/spiceracs/DR1/spice-racs.dr1.corrected.xml


module load rclone

rclone sync -P --transfers $SLURM_NTASKS --checkers $SLURM_NTASKS $data_dir $data_dir_s
rclone sync -P --transfers $SLURM_NTASKS --checkers $SLURM_NTASKS $polcat $polcat_s