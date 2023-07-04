#!/bin/bash -l

#SBATCH --job-name=processSPICE
#SBATCH --export=NONE
#SBATCH --mail-user=alec.thomson@csiro.au
#SBATCH --mail-type=ALL
#SBATCH -e /group/askap/athomson/projects/arrakis/spica/slurmLogs/slurm-%j.err
#SBATCH -o /group/askap/athomson/projects/arrakis/spica/slurmLogs/slurm-%j.out
#SBATCH --ntasks=3

#SBATCH --time=1-00:00:00
#SBATCH --cluster=galaxy
#SBATCH --account=askap


export OMP_NUM_THREADS=1

# Module sorting
conda activate spice
module load singularity
export SINGULARITY_BINDPATH=$(pwd),/group

SPICA=(
    1416+00A
    1351+00A
    1326+00A
    1237+00A
    1302+00A
    1416-06A
    1351-06A
    1326-06A
    1302-06A
    1237-06A
    1418-12A
    1353-12A
    1328-12A
    1303-12A
    1237-12A
    1424-18A
    1357-18A
    1331-18A
    1305-18A
    1239-18A
    1429-25A
    1402-25A
    1335-25A
    1307-25A
    1240-25A
    1212+00A
    1212-06A
    1212-12A
    1213-18A
    1213-25A
)

config=/group/askap/athomson/projects/arrakis/spica/spica_config.txt

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
        data_dir=/group/ja3/athomson/spica

        # Image dirctory
        cd $data_dir

        srun -n 3 --export=ALL processSPICE $field $data_dir/$cal_sbid/RACS_test4_1.05_$field --config $config --savePlots --tt0 $tt0_dir/RACS_test4_1.05_$field.fits --tt1 $tt1_dir/RACS_test4_1.05_$field.fits --use_mpi --skip_cutout --skip_linmos --skip_frion --skip_rmsynth --skip_rmclean

    done
