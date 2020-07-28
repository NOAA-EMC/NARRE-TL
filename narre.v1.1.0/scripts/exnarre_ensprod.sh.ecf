#!/bin/ksh
#####################################################
# produce commom ensemble products (mean, spread and
# probability) of selected variables and write them
# out in grib1 format
#
#  Scirpt: narre_prepare_run_prod.sh.sms
# Purpose: run ensemble product generator to generate 
#          required ensemble products accoriding to variable.tbl 
#  Author: B. Zhou IMSG/EMC/NCEP  10/19/2011 Created
#          Binbin Zhou            02/08/2016 Changed to grib2 
#

set -x

#export NDATE=/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.0/exec/ips/ndate

export XLFRTEOPTS="namelist=old"

narre_cyc=$cyc
grid=$1

#mkdir -m 775 $DATA/ensprod/$grid


wait_time=0
while [ ! -s $DATA_IN/datadone/done4NAMcyc${cyc}grid$grid ] || [ ! -s $DATA_IN/datadone/done4RAPcyc${cyc}grid$grid ] ; do
sleep 10
wait_time=`expr $wait_time + 10`
if [ $wait_time -gt 1800 ] ; then
  echo "Someting wrong for grid" $grid
  exit
fi
done

ls -l $DATA_IN/datadone/done4NAMcyc${cyc}grid$grid $DATA_IN/datadone/done4RAPcyc${cyc}grid$grid

KM=39                                 #vertical pressure levels for all grid#


yy=`echo $PDY | cut -c 1-4`
mm=`echo $PDY | cut -c 5-6`
dd=`echo $PDY | cut -c 7-8`


rundir=$DATA/run.${PDY}/${grid}

# clean directory before job starts
mkdir -p $rundir
cd $rundir
rm -f $rundir/nrre.t*z.m*.* 

 
MODEL_DIR=$DATA_IN/scrub1_${grid}


typeset -Z2 HH
HH=00

MBR[0]=m01
MBR[1]=m02
MBR[2]=m03
MBR[3]=m04
MBR[4]=m05
MBR[5]=m06
MBR[6]=m07
MBR[7]=m08
MBR[8]=m09
MBR[9]=m10
MBR[10]=m11


echo 'narre_cyc=' $narre_cyc
rm -f filename.*

 if   [ $narre_cyc = '00' -o $narre_cyc = '06' -o $narre_cyc = '12' -o $narre_cyc = '18' ] ; then
  times=' 6 12 18 24'
 elif [ $narre_cyc = '01' -o $narre_cyc = '07' -o $narre_cyc = '13' -o $narre_cyc = '19' ] ; then
  times='  7 13 19 25'
 elif [ $narre_cyc = '02' -o $narre_cyc = '08' -o $narre_cyc = '14' -o $narre_cyc = '20' ] ; then
  times='  8 14 20 26'
 elif [ $narre_cyc = '03' -o $narre_cyc = '09' -o $narre_cyc = '15' -o $narre_cyc = '21' ] ; then
  times='  3 9 15 21'
 elif [ $narre_cyc = '04' -o $narre_cyc = '10' -o $narre_cyc = '16' -o $narre_cyc = '22' ] ; then
  times='  4 10 16 22'
 elif [ $narre_cyc = '05' -o $narre_cyc = '11' -o $narre_cyc = '17' -o $narre_cyc = '23' ] ; then
  times='  5 11 17 23'
 fi

NARRE_head=nrre.t${narre_cyc}z
PRCIP_head=prcip.t${narre_cyc}z

if [ $grid = '130' ] ; then
  #fhrs='01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18'
  fhrs='01 02 03 04 05 06 07 08 09 10 11 12'
  ftime=12
else
  fhrs='01 02 03 04 05 06 07 08 09 10 11 12'
  ftime=12
fi

for fhr in $fhrs ; do                                         # for 12 forecast hours

  echo  ${yy} ${mm} ${dd} ${narre_cyc} $fhr $grid $KM $ftime 1 $ftime >> filename.f${fhr}

  MBRS=0                                                            #total number of available members counted for fhr=1

    for time in $times ; do                                            # for 4 NAM members

     NAM_YMDH=`$NDATE -$time ${PDY}${narre_cyc}`
     NAM_dir=`echo $NAM_YMDH | cut -c 1-8`
     NAM_cyc=`echo $NAM_YMDH | cut -c 9-10`

     weight=`echo "scale=2;  1-$time/30" | bc`
     if [ $weight -eq 1.0 ] ; then
        wgt=$weight
     else
        wgt=0$weight
     fi


     HH=`expr $fhr + $time`
     NAM_grib2=$MODEL_DIR/narre.${NAM_dir}/nam.t${NAM_cyc}z.pgrb${grid}.f${HH}.grib2
     echo $NAM_grib2
     if [  -s $NAM_grib2 ] ; then
       ln -sf $NAM_grib2 ${NARRE_head}.${MBR[$MBRS]}.f${fhr}
       #ln -sf ${NAM_grib2}.prcp ${PRCIP_head}.${MBR[$MBRS]}.f${fhr}
       cp ${NAM_grib2}.prcp ${PRCIP_head}.${MBR[$MBRS]}.f${fhr}
       echo "   "$wgt ${NARRE_head}.${MBR[$MBRS]}.f${fhr} '->' $NAM_grib2 >> filename.f${fhr}
       MBRS=`expr $MBRS + 1`
     else
        M1=`expr $MBRS - 1`     #for missing NAM, use previous one as replace
        cp ${NARRE_head}.${MBR[$M1]}.f${fhr} ${NARRE_head}.${MBR[$MBRS]}.f${fhr}
        echo "   "$wgt ${NARRE_head}.${MBR[$MBRS]}.f${fhr} '->' $NAM_grib2 >> filename.f${fhr}
        MBRS=`expr $MBRS + 1`
     fi

    done

  for  time in 0 1 2 3 4 5 6 ; do           # for 7 RAP members, RAP members put last in filename.fxx files
                                            # because in packing cloudbase/ceiling prob, nam's gribfld can not be used (ERROR)
                                            # only rap's cloudbase/ceiling gribfld can be used 

     RAP_YMDH=`$NDATE -$time ${PDY}${narre_cyc}`
     RAP_dir=`echo $RAP_YMDH | cut -c 1-8`
     RAP_cyc=`echo $RAP_YMDH | cut -c 9-10`

     weight=`echo "scale=2;  1-$time/30" | bc`
     if [ $weight -eq 1.0 ] ; then
        wgt=$weight
     else
        wgt=0$weight
     fi

     HH=`expr $fhr + $time`
     RAP_grib2=$MODEL_DIR/narre.${RAP_dir}/rap.t${RAP_cyc}z.pgrb${grid}.f${HH}.grib2
     #if [  -s $RAP_gribr2 ] ; then
     if [  -r $RAP_grib2 ] ; then   #be sure file size exists AND readable as well
      ln -sf $RAP_grib2 ${NARRE_head}.${MBR[$MBRS]}.f${fhr}
      #ln -sf ${RAP_grib2}.prcp  ${PRCIP_head}.${MBR[$MBRS]}.f${fhr}
      cp ${RAP_grib2}.prcp  ${PRCIP_head}.${MBR[$MBRS]}.f${fhr}
      echo  "   "$wgt ${NARRE_head}.${MBR[$MBRS]}.f${fhr} '->' $RAP_grib2 >> filename.f${fhr}
      MBRS=`expr $MBRS + 1`
     fi

  done


done

for mdl in m01 m02 m03 m04 m05 m06 m07 m08 m09 m10 m11 ; do 
 for fhr in $fhrs ; do
  if [ ! -s prcip.t${narre_cyc}z.${mdl}.f${fhr} ] ; then
   rm -f prcip.t${narre_cyc}z.${mdl}.f${fhr}
  fi 
 done 
 echo "${grid} ${mdl} ${narre_cyc}"| $EXECnarre/narre_precip3hr
 for fhr in $fhrs ; do
   cat hr3_prcip.t${narre_cyc}z.${mdl}.f${fhr} >>${PRCIP_head}.${mdl}.f${fhr}
 done
done 


cp $PARMnarre/narre_variable.tbl.grib2.${grid} variable.tbl

for fhr in $fhrs ; do 
  ln -sf filename.f${fhr} filename
  #/nwprod/spa_util/valgrind.v3.11.0/bin/valgrind --leak-check=yes --main-stacksize=20000000 $EXECnarre/narre_ensprod 1>>$pgmout 2>errfile
  $EXECnarre/narre_ensprod 1>>$pgmout 2>errfile
done 

export err=$?;err_chk


if [ $SENDCOM = YES ]; then
 for fhr in $fhrs ; do 
    cp $rundir/nrre.mean.t${narre_cyc}z.f${fhr}     $COMOUT/narre.t${narre_cyc}z.mean.grd${grid}.f${fhr}.grib2
    cp $rundir/nrre.sprd.t${narre_cyc}z.f${fhr}     $COMOUT/narre.t${narre_cyc}z.spread.grd${grid}.f${fhr}.grib2
    cp $rundir/nrre.prob.t${narre_cyc}z.f${fhr}     $COMOUT/narre.t${narre_cyc}z.prob.grd${grid}.f${fhr}.grib2

    if [ $SENDDBN = YES ]; then
      $DBNROOT/bin/dbn_alert MODEL NARRE_GB2 $job $COMOUT/narre.t${cyc}z.mean.grd${grid}.f${fhr}.grib2
      $DBNROOT/bin/dbn_alert MODEL NARRE_GB2 $job $COMOUT/narre.t${cyc}z.spread.grd${grid}.f${fhr}.grib2
      $DBNROOT/bin/dbn_alert MODEL NARRE_GB2 $job $COMOUT/narre.t${cyc}z.prob.grd${grid}.f${fhr}.grib2
    fi
 done
fi

exit
