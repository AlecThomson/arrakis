#!/bin/bash -l

#SBATCH --job-name=processSPICE
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/arrakis/spica/slurmLogs/slurm-%j.err
#SBATCH -o /group/askap/athomson/projects/arrakis/spica/slurmLogs/slurm-%j.out
#SBATCH --ntasks=600
#SBATCH --array=0-30

#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
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

# SPICA=(
#     '1416+00A'
# )

SPICA=(
    '1416+00A'
    '1351+00A'
    '1326+00A'
    '1237+00A'
    '1302+00A'
    '1416-06A'
    '1351-06A'
    '1326-06A'
    '1302-06A'
    '1237-06A'
    '1418-12A'
    '1353-12A'
    '1328-12A'
    '1303-12A'
    '1237-12A'
    '1424-18A'
    '1357-18A'
    '1331-18A'
    '1305-18A'
    '1239-18A'
    '1429-25A'
    '1402-25A'
    '1335-25A'
    '1307-25A'
    '1240-25A'
    '1212+00A'
    '1212-06A'
    '1212-12A'
    '1213-18A'
    '1213-25A'
)

# SPICA=(
#     1213-18A
#     1213-25A
#     1305-18A
#     1307-25A
#     1331-18A
#     1357-18A
#     1402-25A
#     1424-18A
#     1302-06A
#     1351+00A
#     1418-12A
#     1212-06A
#     1237-06A
#     1328-12A
#     1351-06A
#     1416-06A
#     1212-12A
#     1237-12A
#     1303-12A
#     1353-12A
#     1326+00A
#     1416+00A
# )

for field in ${SPICA[*]}
    do
        echo Running pipeline on $field
        cal_sbid=`find_sbid.py $field --cal`
        tt0_dir=/group/askap/athomson/projects/RACS/CI0_mosaic_1.0
        tt1_dir=/group/askap/athomson/projects/RACS/CI1_mosaic_1.0
        # data_dir=/scratch/ja3/athomson/spica
        data_dir=/group/ja3/athomson/spica
        config=/group/askap/athomson/projects/arrakis/spica/spica_config.txt

        # Image dirctory
        cd $data_dir

        # Make a copy of this sbatch file for posterity
        sedstr="s/sbatch/${field}.${SLURM_JOB_ID}\.sbatch/g"
        slurmdir="/group/askap/athomson/projects/arrakis/spica/slurmFiles"
        currentdir="/group/askap/athomson/repos/arrakis/submit"
        sedstr2="s+${currentdir}+${slurmdir}+g"
        # thisfile="/group/askap/athomson/repos/arrakis/submit/processSPICE.sbatch"
        # cp "$thisfile" "$(echo "$thisfile" | sed -e "$sedstr" | sed -e "$sedstr2")"

        # srun --export=ALL processSPICE $field $data_dir/$cal_sbid/RACS_test4_1.05_$field --config $config --savePlots --tt0 $tt0_dir/RACS_test4_1.05_$field.fits --tt1 $tt1_dir/RACS_test4_1.05_$field.fits --use_mpi --skip_cutout
        srun --export=ALL processSPICE $field $data_dir/$cal_sbid/RACS_test4_1.05_$field --config $config --savePlots --tt0 $tt0_dir/RACS_test4_1.05_$field.fits --tt1 $tt1_dir/RACS_test4_1.05_$field.fits --use_mpi --skip_cutout
    done
