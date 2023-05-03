#!/bin/bash -l

#SBATCH --job-name=backup
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/arrakis/spica/slurmLogs/backup_slurm-%j.log
#SBATCH -o /group/askap/athomson/projects/arrakis/spica/slurmLogs/backup_slurm-%j.log
#SBATCH --cluster=zeus
#SBATCH --account=askap
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=2-00:00:00
#SBATCH --partition=copyq

ACA_ALIAS=askap
DATA_DIR=/group/askap/athomson/projects
STAGE_DIR=/scratch/ja3/athomson/staging
BUCKET_NAME=$(basename $DATA_DIR)-bak #-$(date "+%Y-%m-%d")
# Replace underscores with hyphens in bucket name
BUCKET_NAME=${BUCKET_NAME//_/-}

module load rclone


# Make sure the staging directory exists
STAGE_AREA=$STAGE_DIR/$(basename $DATA_DIR)
mkdir $STAGE_AREA

# Get list of files
FILE_LIST=$STAGE_AREA/file_list.tmp
find $DATA_DIR -maxdepth 1 -mindepth 1 -type f > $FILE_LIST
DIR_LIST=$STAGE_AREA/dir_list.tmp
find $DATA_DIR -maxdepth 1 -mindepth 1 -type d > $DIR_LIST

# Copy the data to the staging area
echo "Copying data to staging area"
cat $STAGE_AREA/file_list.tmp | xargs -I {} -P 20 cp -v {} $STAGE_AREA
# tar directories to the staging area
echo "Tarring directories to staging area"

function tar_prog() {
    local DIR=$1
    local STAGE_AREA=$2
    echo $DIR
    tar cf - $DIR -P | pv -s $(du -sb $DIR | awk '{print $1}') > $STAGE_AREA/$(basename $DIR).tar
}
export -f tar_prog
cat $DIR_LIST | xargs -I {} -P 20 bash -c "tar_prog {} $STAGE_AREA"

wait

echo "Uploading data to S3"
# Make sure the bucket exists
rclone mkdir $ACA_ALIAS:$BUCKET_NAME

# Copy the staging area to the bucket
rclone copy $STAGE_AREA $ACA_ALIAS:$BUCKET_NAME --transfers=20 --checkers=20 -P

# Inspect file tree
rclone tree $ACA_ALIAS:$BUCKET_NAME