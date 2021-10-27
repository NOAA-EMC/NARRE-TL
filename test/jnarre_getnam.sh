#!/bin/bash

#PBS -N jnarre_getnam
#PBS -o  /lfs/h2/emc/ptmp/Binbin.Zhou/narre/tmpnwprd1/narre_getnam
#PBS -e  /lfs/h2/emc/ptmp/Binbin.Zhou/narre/tmpnwprd1/narre_getnam
#PBS -q dev
#PBS -S /bin/bash
#PBS -l select=1:ncpus=2:mem=2GB
#PBS -l walltime=01:00:00
#PBS -A NARRE-DEV
#PBS -l debug=true


export narre_ver=v1.2.0
export HOMEnarre=/lfs/h2/emc/vpppg/noscrub/Binbin.Zhou/narre.${narre_ver}

source $HOMEnarre/versions/run.ver
#source $HOMEnarre/modulefiles/$narre_ver 
module load envvar/$envvar_ver
module load PrgEnv-intel/$PrgEnv_intel_ver
module load intel/$intel_ver
module load craype/$craype_ver
module load cray-mpich
module load libjpeg/$libjpeg_ver
module load grib_util/$grib_util_ver
module load wgrib2/$wgrib2_ver
module load cray-pals/$cray_pals_ver

export RUN_ENVIR=dev
export PDY=20210824
export cyc=10

export RAPprod=/lfs/h1/ops/canned/com/rap/v5.1
export COMINnam=/lfs/h1/ops/canned/com/nam/v4.2

$HOMEnarre/jobs/JNARRE_GETNAM $cyc

