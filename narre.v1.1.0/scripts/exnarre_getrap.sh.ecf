#!/bin/ksh
#
#  Script: exnarre_getrap.sh.sms 
#
# purpose:  to get RAP historical files
#           Some addtional aviation products are pre-calculated: CAT, LLWS, 1-hr APCP
#           
#  Author: Binbin Zhou, IMSG/EMC/NCEP
#          10/17/2011    
#          12/10/2012 Binbin Zhou, moved to WCOSS
#          02/08/2016 Binbin Zhou, changed to grib2
#####################################################################################

set -x

#NDATE=${NDATE:-/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.0/exec/ips/ndate}


typeset -Z2 narre_cyc
narre_cyc=$cyc

mkdir -p $DATA_IN/datadone

grid=$1

if [ $grid = '130' ] ; then
  #fhrs='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18'
  fhrs='1 2 3 4 5 6 7 8 9 10 11 12'
else
  fhrs='1 2 3 4 5 6 7 8 9 10 11 12'
fi 

typeset -Z2 HH
typeset -Z2 HH1

echo 'narre_cyc:' $PDY$narre_cyc 

for time in 0 1 2 3 4 5 6 ; do

  RAP_YMDH=`$NDATE -$time ${PDY}${narre_cyc}`
  RAP_dir=`echo $RAP_YMDH | cut -c 1-8`
  RAP_cyc=`echo $RAP_YMDH | cut -c 9-10`
  echo '  Day ' $RAP_dir 'RAP_cyc ' $RAP_cyc

  if [ ${grid} = '130' ] ; then
    RAP_head=$RAPprod/rap.${RAP_dir}/rap.t${RAP_cyc}z.awp${grid}pgrbf
  elif [ ${grid} = '242' ] ; then
    RAP_head=$RAPprod/rap.${RAP_dir}/rap.t${RAP_cyc}z.awp${grid}f
  fi

    scrub1=${DATA_IN}/scrub1_$grid/narre.${RAP_dir}

  if [ ! -d $scrub1 ] ; then
    mkdir -p $scrub1
  fi


  for fhr in $fhrs ; do
   HH=`expr $fhr + $time` 
   RAP_com=${RAP_head}${HH}.grib2
   RAP_scrub1=$scrub1/rap.t${RAP_cyc}z.pgrb${grid}.f${HH}.grib2

   if [ -s $RAP_com ] ; then
     if [ ! -s $RAP_scrub1 ] ; then     
       cd $scrub1
       ln -sf  $RAP_com $RAP_scrub1                                                  #copy rap grib2 files from /com to /scrub2 dir

       #### split 1hr precip data from $RAP_grib1
       HH1=`expr $HH - 1`
       echo "rap.t${RAP_cyc}z.pgrb${grid}.f${HH1}.grib2 rap.t${RAP_cyc}z.pgrb${grid}.f${HH}.grib2 $RAP_cyc $HH1 $HH ${grid} yes ${fhr}"|$EXECnarre/narre_precip > outrap.$fhr

     else
       echo $RAP_scrub1 ' already existing!'
     fi 
    else
     echo $RAP_com ' not existing, skip it!'
    fi

  done   # done for times 

done     # done for fhr

    echo 'Done for RAP cyc ' $cyc ' grid' $grid  >> $DATA_IN/datadone/done4RAPcyc${cyc}grid$grid

    ls -l $DATA_IN/datadone/done4RAPcyc${cyc}grid$grid

