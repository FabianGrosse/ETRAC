#!/bin/bash
#
# Project identification
#
# SBATCH --account=JOBACC_TMP #no need, this is setup in ~/.bashrc
#
# MAX. JOB LENGTH
#
#SBATCH --time=JOBTIME_TMP
#
# NODES, CORES/NODE
#
#SBATCH --ntasks=NTASKS_TMP
#SBATCH --cpus-per-task=1
#
#SBATCH --mem-per-cpu=MEM_CPU_TMP # memory; default unit is megabytes

# ENVIRONMENT VARIABLES
#
#SBATCH --export=ALL
#
# DO NOT RESTART AUTOMATICALLY
#
#SBATCH --no-requeue
#
# JOB NAME
#
#SBATCH --job-name=JOBNAME_TMP
#
# EMAIL JOB RESULTS
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user="USER_EMAIL_TMP"
#
# LOG FILE NAMES
#
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# LOAD NETCDF LIBRARY
module load netcdf-fortran/4.4.4

# SEND JOB
echo ' '
echo '======================================================'
echo '  ETRAC start running at '`date`
echo '------------------------------------------------------'
# tell wrapper script that job started
# => required for time tracking within wrapper script
echo "job started" > jobStart
srun ./ETRAC
# tell wrapper job that job ended
echo "job ended" > jobEnd
echo '------------------------------------------------------'
echo ' ETRAC finished at      '`date`
echo '======================================================'
echo ' '

exit
