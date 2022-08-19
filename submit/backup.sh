#!/bin/bash -l

#SBATCH --job-name=backup
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/spiceracs/spica/slurmLogs/backup_slurm-%j.log
#SBATCH -o /group/askap/athomson/projects/spiceracs/spica/slurmLogs/backup_slurm-%j.log
#SBATCH --cluster=zeus
#SBATCH --account=askap
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=2-00:00:00
#SBATCH --partition=copyq

ACA_ALIAS=acacia-spiceracs
DATA_DIR=/group/ja3/athomson/full_spica
BUCKET_NAME=athomson-fullspica

module load rclone
# rclone sync -P -L $data_dir acacia-spiceracs:$BUCKET_NAME --transfers=20 --checkers=20
# tar -zcvf - $DATA_DIR | rclone rcat acacia-spiceracs:$BUCKET_NAME/data_date "+%Y-%m-%d_%H:%M:%S".tar.gz -v
# rclone tree acacia-spiceracs:$BUCKET_NAME

# Make sure the bucket exists
rclone mkdir $ACA_ALIAS:$BUCKET_NAME

# Copy files to the bucket
for f in $(find $DATA_DIR -maxdepth 1 -type f)
do
    rclone copy $f $ACA_ALIAS:$BUCKET_NAME/ -P &
done

# Copy directories to the bucket and tar them on the way
for f in $(find $DATA_DIR -maxdepth 1 -type d)
do
    tar -zcvf - $f | rclone rcat $ACA_ALIAS:$BUCKET_NAME/$(basename $f).tar.gz -v &
done

wait

rclone tree acacia-spiceracs:$BUCKET_NAME