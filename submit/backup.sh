#!/bin/bash -l

#SBATCH --job-name=backup
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/backup_slurm-%j.err
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/backup_slurm-%j.out

#SBATCH --cluster=zeus
#SBATCH --account=askap
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=2-00:00:00
#SBATCH --partition=copyq


data_dir=/group/ja3/athomson/spica
tar_name=spica.tar.gz

module load rclone
tar -cf - $data_dir | pigz | rclone rcat --s3-chunk-size 600M acacia-spiceracs:spica/$tar_name  --progress
rclone tree acacia-spiceracs: