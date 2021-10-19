#!/bin/bash

#PBS -N jnarre_getnam
#PBS -o  /lfs/h2/emc/ptmp/Binbin.Zhou/narre/tmpnwprd1/narre_getnam
#PBS -e  /lfs/h2/emc/ptmp/Binbin.Zhou/narre/tmpnwprd1/narre_getnam
#PBS -q dev
#PBS -l select=1:ncpus=4
#PBS -l walltime=01:30:00
#PBS -A NARRE-DEV


export narre_ver=v1.2.0
export HOMEnarre=/lfs/h2/emc/vpppg/noscrub/Binbin.Zhou/narre.${narre_ver}

source $HOMEnarre/versions/run.ver
source $HOMEnarre/modulefiles/$narre_ver 


export RUN_ENVIR=dev
export PDY=20210824
export cyc=05

export RAPprod=/lfs/h1/ops/canned/com/rap/v5.1
export COMINnam=/lfs/h1/ops/canned/com/nam/v4.2

$HOMEnarre/jobs/JNARRE_GETNAM 
