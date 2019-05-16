#!/bin/bash
# =====================================================================
# script for running a ETRAC calculation as a series of subsequent jobs
# by Fabian Grosse (fabian.grosse[at]dal.ca)
# =====================================================================
# relevant functions:
# - automatic continuation if job was cancelled due to cluster failure
#
# - email notifications are send to user if:
#   1) a (sub-)job fails before entering time loop => job abortion
#   2) a simulation blows up                       => job abortion
#   3) a (sub-)job failed due to cluster failure   => auto-continuation
#        or time-out
# =====================================================================
# HOW TO USE?
# 1) define set-up in script header
# 2) execute script by: nohup ./run_ETRAC_CHAIN.sh &
# =====================================================================
# last changed: May 16th, 2019
# =====================================================================
# USER-DEFINED SET-UP
# =====================================================================

# debug this script? yes (1) or no (0 = default)?
# - debug=1: only tests this script; script log is prompted to terminal, no job is sumbitted
# - debug=0: script log is written to log file
debug=0

# define calculation ID (used for names of output files)
#runID=roms854-NGoMex-TBNT
runID=roms854-MCH-ETRAC

# set directories and names of input, output and temporary files
# NOTE: directories require '/' at end of name
#
# path and file name of merged yearly ROMS output files
# for running multiple years, use YEAR_TMP placeholder in file name
bulkPath=/project/def-kfennel/grosse/ETRAC_data/
bulkDummy=tbnt_mch_bio_tbnt_YEAR_TMP.nc
# path and file name of initialization files (i.e. distributions of relative tracer fractions)
#  - for running multiple years, use INIYEAR_TMP placeholder in file name
#  - init file name must be of structure "${runID}_INIYEAR_TMP_relative_fractions.nc"
#    when running a sequence of multiple years
initPath=/project/def-kfennel/grosse/ETRAC_data/
initDummy=${runID}_INIYEAR_TMP_relative_fractions.nc
# output path
#  - initPath and outputPath must match in case of a sequence of jobs
outputPath=/project/def-kfennel/grosse/ETRAC_data/
# temporary path used only during the job
tmpPath=/scratch/grosse/ETRAC_${runID}_TMP/

# project and work paths
#  - the project path contains all scripts (incl. this very one), template files, the 'Build' folder
#    with the compiled executable and the folders 'model_setup' and 'etrac_setup'
#  - this script needs to be executed from the project directory
#  - copies of the jobs' run files and the log file are stored in the work path
projPath=`pwd`
wrkPath=${projPath}/wrk

# directories with model and tracing setup files used by ETRAC
# NOTE: directories require '/' at end of name
modelPath=/home/grosse/ETRAC/setup_files/model_setup/
etracPath=/home/grosse/ETRAC/setup_files/etrac_setup/

# set template file used for ETRAC setup
setBase=etrac_set_BASE.nml

# define file name for log file => must be identical to the name defined in ETRAC/software/src/etrac_common.F90
logFile=etrac_logfile.dat

# define if job is a new job, i.e., starting from initialization (0) or not (1)
# set to 1 if you start an entirely new job, set to 0 if you have an initialization file available
# or if you have to manually continue a failed job
newJob=0

# submit job or only prepare setup? (only useful if single job)
submitJob=1

# keep log files of individual jobs? yes (1) or no (0)
keepLog=1

# keep setup and slurm script files of individual jobs? yes (1) or no (0)
keepRunFiles=0

# define time step and sub-time step (in minutes), and time step index for first step of first job
#  - as ETRAC applies an automatic recursive time step subdivision (if needed), SUB_STEP should be
#    identical to TIME_STEP; if you want to apply a fixed smaller sub step, SUB_STEP must be a proper
#    divider of TIME_STEP (note that automatic time step subdivision may still occur)
#  - set 'FINAL_STEP=0' to run until end of each year
TIME_STEP=1440
SUB_STEP=1440
FIRST_STEP=50
FINAL_STEP=51

# if 1st year's calculation is not started on very first time step, define time step offset
# (i.e., actual starting day) to ensure correct reading of bulk variables
OFFSET_STEP=49

# provide start and end year
FIRSTYEAR=2001
LASTYEAR=2001

# define SLURM job settings
# (account name, memory per CPU, job name, ntasks, email of user)
JOB_ACC='def-kfennel'
MEM_CPU='7000M'
JOB_NAME=ETRAC_${runID}_CHAIN
NTASKS=1
USER_EMAIL='DUMMY.USER@EMAIL.COM'
# time limit: days, hours, minutes and seconds; individual 2-digit strings
JOB_DD='02' # days
JOB_HH='23' # hours
JOB_MM='59' # minutes
JOB_SS='55' # seconds

# ===================================================================
# END OF USER-DEFINED SET-UP
# ===================================================================
# AUTOMATIC JOB SETUP AND SUBMISSION
# ===================================================================

# filename of this script's log
scriptLog=CHAIN_LOG_${runID}.log

# define default job time string and value in seconds
JOB_TIME_STR=${JOB_DD}-${JOB_HH}:${JOB_MM}:${JOB_SS}
JOB_TIME_SEC=$((10#${JOB_DD}*86400 + 10#${JOB_HH}*3600 + 10#${JOB_MM}*60 + 10#${JOB_SS}))

# calculate overall number of jobs
if [ ${FINAL_STEP} -le 0 ]; then
  let NJOBS=LASTYEAR-FIRSTYEAR+1
else 
  NJOBS=1
fi

# get digits of maximum job number
if [ ${NJOBS} -lt 10 ]; then
  JOBDIG=1
else
  if [ ${NJOBS} -lt 100 ]; then
    JOBDIG=2
  else
    if [ ${NJOBS} -lt 1000 ]; then
      JOBDIG=3
    else
      echo "Maximum number of consecutive jobs for one simulation is 999."
      exit
    fi
  fi
fi

# create output directory
mkdir -p ${outputPath}

# create work directory
mkdir -p ${wrkPath}

# start logging 
rm -f ${scriptLog}
head="=======================\n CHAIN JOB INFORMATION\n======================="
msg1=">> ETRAC setup <<\n runID = ${runID}\n newJob = ${newJob}\n bulk path = ${bulkPath}\n bulk dummy = ${bulkDummy}\n init path = ${initPath}\n init dummy = ${initDummy}\n output path = ${outputPath}\n temporary path = ${tmpPath}\n log file = ${logFile}\n"
msg2=">> chain job setup <<"
if [ ${debug} -eq 1 ]; then
  echo -e ${head}
  echo -e ${msg1}
  echo -e ${msg2}
else
  echo -e ${head} > ${scriptLog}
  echo -e ${msg1} >> ${scriptLog}
  echo -e ${msg2} >> ${scriptLog}
fi

# update non-time-varying placeholders in etrac_set_BASE.nml and store as new file
awk '{gsub("RUNID_TMP", "'${runID}'"); \
      gsub("TEMPPATH_TMP", "'${tmpPath}'"); \
      gsub("BULKPATH_TMP", "'${bulkPath}'"); \
      gsub("BULKFILE_TMP", "'${bulkDummy}'"); \
      gsub("INITPATH_TMP", "'${initPath}'"); \
      gsub("INITFILE_TMP", "'${initDummy}'"); \
      gsub("OUTPUTPATH_TMP", "'${outputPath}'"); \
      gsub("DT_TMP", "'${TIME_STEP}'"); \
      gsub("DT_SUB_TMP", "'${SUB_STEP}'"); \
      print}' ${setBase} > etrac_set_${runID}_BASE.nml

# update non-time-varying placeholders in SLURM run script (run_ROMS_SINGLE_BASE.sh) and store as new file
awk '{gsub("JOBACC_TMP", "'${JOB_ACC}'"); \
      gsub("NTASKS_TMP", "'${NTASKS}'"); \
      gsub("MEM_CPU_TMP", "'${MEM_CPU}'"); \
      gsub("JOBTIME_TMP", "'${JOB_TIME_STR}'"); \
      gsub("JOBNAME_TMP", "'${JOB_NAME}'"); \
      gsub("USER_EMAIL_TMP", "'${USER_EMAIL}'"); \
      print}' ETRAC_SINGLE_BASE.slurm > ETRAC_${runID}_BASE.slurm


# loop over jobs: placeholders in template files are updated and jobs are submitted successively
IJOB=1
jobOK=1
RSTcount=0

rm -rf ${tmpPath} ${wrkPath}${logFile}

while [ ${IJOB} -le ${NJOBS} ]; do
  
  # set start and end step of job
  if [ ${jobOK} -eq 1 ]; then
    # previous job succeeded => start calculation for next year
    let YEAR=FIRSTYEAR+IJOB-1
    if [ ${FINAL_STEP} -gt 0 ]; then
      NSTEPS=${FINAL_STEP}
    else
      # leap year check
      if [ `expr ${YEAR} % 400` -eq 0 ]; then
        let NSTEPS=366*1440/TIME_STEP
      else
        if [ `expr ${YEAR} % 4` -eq 0 ] && [ `expr ${YEAR} % 100` -ne 0 ]; then
          let NSTEPS=366*1440/TIME_STEP
        else
          let NSTEPS=365*1440/TIME_STEP
        fi
      fi
    fi
    if [ ${debug} -eq 1 ]; then
      echo $IJOB $YEAR $NSTEPS
    fi
    if [ ${newJob} -eq 1 ]; then
      # new job: 1st job starts from arbitrary initial conditions; subsequent jobs use results of previous job
      if [ ${IJOB} -eq 1 ]; then
        isWarm=0
        startStep=${FIRST_STEP}
	INIYEAR=${YEAR}
      else
        isWarm=1
        startStep=1
	let INIYEAR=YEAR-1
      fi
      contWrt=0
      offset=0
    else
      # continued job using existing results for initialisation
      isWarm=1
      if [ ${IJOB} -eq 1 ]; then
        startStep=${FIRST_STEP}
	if [ ${startStep} -eq 1 ]; then
          let INIYEAR=YEAR-1
	  contWrt=0
        else
          INIYEAR=${YEAR}
          contWrt=1
        fi
      else
        contWrt=0
        startStep=1
        let INIYEAR=YEAR-1
      fi
      if [ ${YEAR} -eq ${FIRSTYEAR} ] && [ ${OFFSET_STEP} -ne 0 ]; then
	offset=${OFFSET_STEP}
      else
	offset=0
      fi
    fi
  else
    # previous job failed => continue calculation
    isWarm=1
    contWrt=1
    INIYEAR=${YEAR}
    startStep=${lastStep}
    if [ ${YEAR} -eq ${FIRSTYEAR} ] && [ ${OFFSET_STEP} -ne 0 ]; then
      offset=${OFFSET_STEP}
    else
      offset=0
    fi
  fi
  endStep=${NSTEPS}
  JOB_TIME=${JOB_TIME_STR}

  # double check if time range is reasonable
  if [ ${startStep} -gt ${endStep} ]; then
    # start step > end step
    # send notification email to user and abort script => no further job submission
    if [ ${debug} -eq 1 ]; then
      echo -e "${runID}_${IJOB}-${RSTcount} not submitted\nInvalid time range (start > end). Check your setup and chain job script!"
    else
      echo -e "Invalid time range (start > end). Check your setup and chain job script!" | mail -s "${runID}_${IJOB}-${RSTcount} not submitted" ${USER_EMAIL}
    fi
    exit
  fi 

  # set file name of ocean.in and run script files
  if [ ${JOBDIG} -eq 1 ]; then
    setFile=etrac_set_${runID}_${IJOB}_${RSTcount}.nml
    runFile=ETRAC_${runID}_${IJOB}_${RSTcount}.slurm
  else
    if [ ${JOBDIG} -eq 2 ]; then
      if [ ${IJOB} -lt 10 ]; then
        setFile=etrac_set_${runID}_0${IJOB}_${RSTcount}.nml
        runFile=ETRAC_${runID}_0${IJOB}_${RSTcount}.slurm
      else
        setFile=etrac_set_${runID}_${IJOB}_${RSTcount}.nml
        runFile=ETRAC_${runID}_${IJOB}_${RSTcount}.slurm
      fi
    else
      if [ ${IJOB} -lt 10 ]; then
        setFile=etrac_set_${runID}_00${IJOB}_${RSTcount}.nml
        runFile=ETRAC_${runID}_00${IJOB}_${RSTcount}.slurm
      else
        if [ ${IJOB} -lt 100 ]; then
          setFile=etrac_set_${runID}_0${IJOB}_${RSTcount}.nml
          runFile=ETRAC_${runID}_0${IJOB}_${RSTcount}.slurm
        else
          setFile=etrac_set_${runID}_${IJOB}_${RSTcount}.nml
          runFile=ETRAC_${runID}_${IJOB}_${RSTcount}.slurm
        fi
      fi
    fi
  fi
  # write job time limit to log file
  msg=" job: ${IJOB}-${RSTcount}\n isWarmStart = ${isWarm}\n continueWrite = ${contWrt}\n run year = ${YEAR}\n init year = ${INIYEAR}\n steps = ${startStep}-${endStep}\n setFile = ${setFile}\n run script = ${runFile}"
  if [ ${debug} -eq 1 ]; then
    echo -e ${msg}
  else
    echo -e ${msg} >> ${scriptLog}
  fi
  # update warm start switch, continued writing switch, years and start and end steps in etrac_set_[runID]_BASE.nml and store as new file for the job 
  awk '{gsub("ISWARM_TMP", "'${isWarm}'"); \
        gsub("CONTINUE_TMP", "'${contWrt}'"); \
        gsub("INIYEAR_TMP", "'${INIYEAR}'"); \
        gsub("YEAR_TMP", "'${YEAR}'"); \
        gsub("START_TMP", "'${startStep}'"); \
        gsub("END_TMP", "'${endStep}'"); \
        gsub("OFFSET_TMP", "'${offset}'"); \
        print}' etrac_set_${runID}_BASE.nml > ${setFile}

  if [ ${debug} -eq 1 ]; then
    # run "fake" job
    let slurmJobID=10*IJOB+RSTcount
    echo -e " slurmJobID = ${slurmJobID}"
    echo "Sleep 5s to mimic job ..."
    startTime=`date +%s`
    sleep 5
    let waitTime="$(date +%s)"-startTime
    echo "Time waited: "${waitTime}"sec"
    cp -f ETRAC_${runID}_BASE.slurm ${runFile}
    if [ ${keepRunFiles} -eq 1 ]; then
      cp -f ${runFile} ${wrkPath}/ETRAC-${slurmJobID}.slurm
      cp -f ${setFile} ${wrkPath}/etrac_set_${slurmJobID}.nml
    fi
  else
    # create temporary work directory and copy relevant files
    mkdir -p $tmpPath
    cd $tmpPath
    cp -f $projPath/ETRAC_${runID}_BASE.slurm ${runFile}
    cp -f $projPath/${setFile} .
    ln -s ${setFile} etrac_set.nml
    cp -f $projPath/Build/ETRAC .
    cp -rf $projPath/$modelPath .
    cp -rf $projPath/$etracPath .
    cp -f $projPath/${scriptLog} .
    # submit job
    if [ ${submitJob} -eq 1 ]; then
      sbatch ${runFile}
    else
      echo "Run files created and copied to: "${tmpPath}". Ready for manual submission."
      exit
    fi
    # sleep until job has started (i.e., until "jobStart" file exists)
    running=0
    while [ ${running} -eq 0 ]; do
      if [ -e jobStart ]; then
        running=1
      else
        sleep 300
      fi
    done
    # check for early job error
    if [ -e jobEnd ]; then
       cp -f etrac_logfile.dat ${wrkPath}${logFile}
       echo -e "ETRAC calculation using ${runFile} and ${setFile} failed prematurely.\nCheck your log file: ${logFile}!" >> ${scriptLog}
       echo -e "Simulation using ${runFile} and ${setFile} failed prematurely.\nCheck your log file: ${logFile}!" | mail -s "${SIMID}_${IJOB}-${RSTcount} failed" ${USER_EMAIL}
       exit
    fi
    # get job ID, start time and remaining job time
    gotInfo=0
    while [ ${gotInfo} -eq 0 ]; do
      squeueOut=`squeue --name=${JOB_NAME} -u ${USER} -o %A,%L,%S,%u | tail -n1`
      gotInfo="$(echo "${squeueOut}" | grep ${USER} | wc -l)"
      if [ ${gotInfo} -eq 1 ]; then
        slurmJobID=`echo "${squeueOut}" | cut -d ',' -f 1`
        leftTimeStr=`echo "${squeueOut}" | cut -d ',' -f 2`
        startTimeStr=`echo "${squeueOut}" | cut -d ',' -f 3`
      else
        sleep 300
      fi
    done
    isDays="$(echo "${leftTimeStr}" | grep - | wc -l)"
    if [ ${isDays} -eq 1 ]; then
      leftDays="$(echo "${leftTimeStr}" | cut -d '-' -f 1)"
      leftTimeStr="$(echo "${leftTimeStr}" | cut -d '-' -f 2)"
    else
      leftDays=0
    fi
    leftHours="$(echo "${leftTimeStr}" | cut -d ':' -f 1)"
    leftMins="$(echo "${leftTimeStr}" | cut -d ':' -f 2)"
    leftSecs="$(echo "${leftTimeStr}" | cut -d ':' -f 3)"
    
    # write start time to file
    echo -e "Job ${slurmJobID} started: ${startTimeStr}" >> ${scriptLog}

    # JOB STARTED => set run time counter to remaining job time
    runTime=$((10#${leftDays}*86400+10#${leftHours}*3600+10#${leftMins}*60+10#${leftSecs}))
    # sleep until job has finished (i.e. until "jobEnd" file exists or "time limit + 2min" is exceeded)
    while [ ${running} -eq 1 ]; do
      if [ -e jobEnd ] || [ ${runTime} -lt -120 ]; then
        running=0
        echo -e "Job ${slurmJobID} ended:   `date`" >> ${scriptLog}
      else
        sleep 300
        let runTime=runTime-300
      fi
    done
    # job complete: copy files to wrk directory
    cd $wrkPath
    mv -f $tmpPath/${logFile} .
    if [ ${keepLog} -eq 1 ]; then
      cp -f  ${logFile} ${logFile}_${slurmJobID}
    fi
    mv -f $tmpPath/${JOB_NAME}-${slurmJobID}* .
    mv -f $tmpPath${runFile} .
    mv -f $tmpPath${setFile} .
    if [ ${keepRunFiles} -eq 1 ]; then
      cp -f ${runFile} ETRAC-${slurmJobID}.slurm
      cp -f ${setFile} etrac_set_${slurmJobID}.nml
    fi
    cd $projPath
    mv -f $tmpPath/${scriptLog} .
    rm -rf ${setFile} $tmpPath
  fi
  
  # after job finished: check log file
  if [ -e ${wrkPath}/${logFile} ]; then
    exist=1 # all good, file exists => just du something, otherwise syntax error
  else
    # log file does not exist
    # send notification email to user and abort script => no further job submission
    if [ ${debug} -eq 1 ]; then
      echo -e "${runID}_${IJOB}-${RSTcount} (slurm ID: ${slurmJobID}) failed\nSimulation using ${runFile} and ${setFile} failed without writing a log file.\nCheck your setup!"
    else
      echo -e "Simulation using ${runFile} and ${setFile} failed without writing a log file.\nCheck your setup!" | mail -s "${runID}_${IJOB}-${RSTcount} (slurm ID: ${slurmJobID}) failed" ${USER_EMAIL}
    fi
    exit
  fi

  if [ "$(grep normal\ end ${wrkPath}/${logFile} | wc -l)" -eq 1 ]; then
    # job succeeded
    jobOK=1
    let IJOB=IJOB+1
    RSTcount=0
    rm -f ${wrkPath}/${runFile} ${wrkPath}/${setFile} ${wrkPath}/${logFile}
    if [ ${IJOB} -gt ${NJOBS} ]; then
      exit
    fi
  else
    # job failed => check if error occured
    if [ "$(grep \ STOP\  ${wrkPath}/${logFile} | wc -l)" -eq 1 ]; then
      # job failed due to unexpected error during execution
      # send notification email to user and abort script => no further job submission
      if [ ${debug} -eq 1 ]; then
        echo -e "${runID}_${IJOB}-${RSTcount} (slurm ID: ${slurmJobID}) failed\nSimulation using ${runFile} and ${setFile} failed due to early program abortion.\nCheck your setup and log file: ${logFile}!"
      else
        echo -e "Simulation using ${runFile} and ${setFile} failed due to early program abortion.\nCheck your setup and log file: ${logFile}!" | mail -s "${runID}_${IJOB}-${RSTcount} (slurm ID: ${slurmJobID}) failed" ${USER_EMAIL}
      fi
      exit
    else
      # job failed due to time limit exceedance
      # check if time loop was entered
      if [ "$(grep MAIN\ STEP ${wrkPath}/${logFile} | wc -l)" -eq 0 ]; then
	# time limit exceeded before entering time loop
	# send notification email to user and abort script => no further job submission
	if [ ${debug} -eq 1 ]; then
          echo -e "${runID}_${IJOB}-${RSTcount} (slurm ID: ${slurmJobID}) failed\nJob time limit exceeded before entering time loop.\nIncrease your time limit!"
        else
          echo -e "Job time limit exceeded before entering time loop.\nIncrease your time limit!" | mail -s "${runID}_${IJOB}-${RSTcount} (slurm ID: ${slurmJobID}) failed" ${USER_EMAIL}
        fi
        exit
      fi
      # time limit exceeded during calculation
      jobOK=0
      # send notification email to user and continue job
      if [ ${debug} -eq 1 ]; then
        echo -e "${runID}_${IJOB}-${RSTcount} (slurm ID: ${slurmJobID}) failed\nSimulation using ${runFile} and ${setFile} failed before anticipated end - time limit exceeded.\nAutomatic re-submission."
      else
        echo -e "Simulation using ${runFile} and ${setFile} failed before anticipated end - time limit exceeded.\nAutomatic re-submission." | mail -s "${runID}_${IJOB}-${RSTcount} (slurm ID: ${slurmJobID}) failed" ${USER_EMAIL}
      fi
      # get time step during which execution stopped
      let lastStep="$(grep MAIN\ STEP: ${wrkPath}/${logFile} | tail -n1 | cut -d ':' -f 2 | sed -e 's/^[[:space:]]*//')"+startStep-1
      let RSTcount=RSTcount+1
      rm -f ${wrkPath}/${runFile} ${wrkPath}/${setFile} ${wrkPath}/${logFile}
      if [ ${debug} -eq 1 ]; then
        echo "lastStep = ${lastStep}"
      fi
    fi
  fi
  
  if [ ${debug} -eq 1 ]; then
    echo -e "jobOK = ${jobOK}\n"
  fi

done

rm -f etrac_set_${runID}_BASE.nml ETRAC_${runID}_BASE.slurm

exit
