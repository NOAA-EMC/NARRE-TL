#!/bin/ksh
#
#  Script:  get_nam_on_rap.sh
#
# purpose:  to get NAM files and copygb them onto RAP CONUS/130 or Alaska/242 grid
#           Some addtional aviation products are pre-calculated: CAT, LLWS, 1-hr APCP
#           
#  Author: Binbin Zhou, IMSG/EMC/NCEP
#          10/17/2011    
#          12/10/2012 Binbin Zhou, moved to WCOSS
#          02/08/2016 Binbin Zhou, changed to grib2
#          12/02/2018 Binbin Zhou, moved to Dell
#####################################################################################

set -x 

#module load ics


#NDATE=/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.0/exec/ips/ndate
typeset -Z2 narre_cyc
narre_cyc=$cyc
grd=$1

if [ $grd = '130' ] ; then
  #fhrs='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18'
  fhrs='1 2 3 4 5 6 7 8 9 10 11 12'
else
  fhrs='1 2 3 4 5 6 7 8 9 10 11 12'
fi 

mkdir -p $DATA_IN/datadone

typeset -Z2 HH
typeset -Z2 HH1
if [ ${grd} = '130' ] ; then
 #grid="130 3 451 337 16281 -126138  8 -95000 13545 13545 0 64   0 25000 25000 0 0"                 #to 13km RAP grid
 grid="30 6 0 0 0 0 0 0 451 337 16281000 233862000 8 25000000 265000000 13545000 13545000 0 64 25000000 25000000 0 0" #grid 130 RAP grib2
 gridwesley="lambert:265.0:25.0:25.0:25.0 233.862:451:13545.0 16.281:337:13545.0"
else
 grid="20 6 0 0 0 0 0 0 553 425 30000000 187000000 8 60000000 225000000 11250000 11250000 0 64"                 #grid 242 RAP grib2
 gridwesley="nps:225:60 187:553:11250 30:425:11250"
fi

 #times is previous hours (ie forecaste hours of NAM).   
 if   [ $narre_cyc = '00' -o $narre_cyc = '06' -o $narre_cyc = '12' -o $narre_cyc = '18' ] ; then  
  times=' 6 12 18 24'
 elif [ $narre_cyc = '03' -o $narre_cyc = '09' -o $narre_cyc = '15' -o $narre_cyc = '21' ] ; then
  times='  3 9 15 21'
 elif [ $narre_cyc = '04' -o $narre_cyc = '10' -o $narre_cyc = '16' -o $narre_cyc = '22' ] ; then
  times='  4 10 16 22' 
 elif [ $narre_cyc = '05' -o $narre_cyc = '11' -o $narre_cyc = '17' -o $narre_cyc = '23' ] ; then
  times='  5 11 17 23'
 elif [ $narre_cyc = '01' -o $narre_cyc = '07' -o $narre_cyc = '13' -o $narre_cyc = '19' ] ; then
  times='  7 13 19 25' 
 elif [ $narre_cyc = '02' -o $narre_cyc = '08' -o $narre_cyc = '14' -o $narre_cyc = '20' ] ; then
  times='  8 14 20 26'
 fi

 echo 'narre_cyc:' $PDY $narre_cyc 'times:' $times

for time in $times ; do

  NAM_YMDH=`$NDATE -$time ${PDY}${narre_cyc}`
  NAM_dir=`echo $NAM_YMDH | cut -c 1-8`
  NAM_cyc=`echo $NAM_YMDH | cut -c 9-10`

  if [ ${grd} = '130' ] ; then
    NAM_com_head=$COMINnam/nam.${NAM_dir}/nam.t${NAM_cyc}z.awphys                    
    NAM_com_awip32=$COMINnam/nam.${NAM_dir}/nam.t${NAM_cyc}z.awip32                  #to get cloud base/top data
  elif [ ${grd} = '242' ]; then
    NAM_com_head=$COMINnam/nam.${NAM_dir}/nam.t${NAM_cyc}z.awip32
  fi

   scrub1=${DATA_IN}/scrub1_$grd/narre.${NAM_dir}

  if [ ! -d $scrub1 ] ; then
    mkdir -p $scrub1
  fi

  cd $scrub1


  for fhr in $fhrs  ; do    # fhr is narre's forecast hours

   HH=`expr $fhr + $time` 
   if [ ${grd} = '130' ] ; then
     NAM_com_awphys=${NAM_com_head}${HH}.tm00.grib2 
   else
     NAM_com_awphys=$COMINnam/nam.${NAM_dir}/nam.t${NAM_cyc}z.awip32${HH}.tm00.grib2
     #$CNVGRIB -g12 -p40 ${NAM_com_head}${HH}.tm00 $NAM_com_awphys         #NAM 221 grid grib2 files only have f00 f03 f06 .. files, so get grib2 files from hourly grib1  
   fi 

   if [ -s $NAM_com_awphys ] ; then 
     

       if [ ${grd} = '130' ] ; then         #Note: for NAM-218: no cloud base/top, so retrieve cloud from awip32 and append it to awphys
          rm -f awip32_cloud.f${HH}.grib2   
          $WGRIB2 -match "HGT:cloud top" $NAM_com_awip32${HH}.tm00.grib2 |$WGRIB2 -i $NAM_com_awip32${HH}.tm00.grib2 -grib cloud
          cat cloud > awip32_cloud.f${HH}.grib2
          $WGRIB2 -match "HGT:cloud base" $NAM_com_awip32${HH}.tm00.grib2 |$WGRIB2 -i $NAM_com_awip32${HH}.tm00.grib2 -grib cloud
          cat cloud >> awip32_cloud.f${HH}.grib2
       fi

       NAM_narre=$scrub1/nam.t${NAM_cyc}z.pgrb${grd}.f${HH}.grib2
  
       if [ ! -s $NAM_narre ] ; then
      
          $WGRIB2 $NAM_com_awphys -new_grid_winds earth -new_grid $gridwesley $NAM_narre 


         if [ ${grd} = '130' ] ; then       #only for NAM-218, there is no cloud base/top, so append it, but first convert to 130 RAP grid         
            $WGRIB2 awip32_cloud.f${HH}.grib2 -set_grib_type same -new_grid_winds grid -new_grid $gridwesley awip32_cloud${grd}.f${HH}.grib2   
            cat   awip32_cloud${grd}.f${HH}.grib2  >> $NAM_narre
         fi

         ####Note:  split 1hr precip data from $NAM_grib_narre
         HH1=`expr $HH - 1` 
         echo "nam.t${NAM_cyc}z.pgrb${grd}.f${HH1}.grib2 nam.t${NAM_cyc}z.pgrb${grd}.f${HH}.grib2 $NAM_cyc $HH1 $HH ${grd} yes ${fhr}" |$EXECnarre/narre_precip

       else
         echo $NAM_grib_narre ' already existing!'
       fi

    else
       echo $NAM_grib ' not existing, skip it!'
    fi
 
  done   

done        

echo 'finish for NAM cyc ' $cyc ' grid' $grd  >> $DATA_IN/datadone/done4NAMcyc${cyc}grid$grd 
ls -l $DATA_IN/datadone/done4NAMcyc${cyc}grid$grd

if [ $grd = 130 ]; then
  ecflow_client --event release_ensprod
fi

exit

