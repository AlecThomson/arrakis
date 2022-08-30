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
DATA_DIR=/group/ja3/athomson/spica
STAGE_DIR=/scratch/ja3/athomson/staging
BUCKET_NAME=$(basename $DATA_DIR)-bak-$(date "+%Y-%m-%d")
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
for DIR in $(cat $DIR_LIST)
do
    # ((i=i%N)); ((i++==0)) && wait
    # tar -cfh $STAGE_AREA/$(basename $DIR).tar $DIR &
    echo $DIR
    tar cf - $DIR -P | pv -s $(du -sb $DIR | awk '{print $1}') > $STAGE_AREA/$(basename $DIR).tar
done

wait

echo "Uploading data to S3"
# Make sure the bucket exists
echo rclone mkdir $ACA_ALIAS:$BUCKET_NAME

# # Copy the staging area to the bucket
rclone copy $STAGE_AREA $ACA_ALIAS:$BUCKET_NAME --transfers=20 --checkers=20

rclone tree $ACA_ALIAS:$BUCKET_NAME