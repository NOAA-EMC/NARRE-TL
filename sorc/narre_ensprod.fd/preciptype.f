ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Subroutine preciptype: Compute precipitation type rain/freezing rain/snow
c  based on Jun Du's old version
c
c  Author: Binbin Zhou, Aug. 4, 2005 
c  Modification history:
c   01/12/2006: Geoff Manikin: Add Geoff Manikin's algorithm to determine dominant precip_type 
c   04/10/2013: Binbin Z. Modified for grib2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   	subroutine  preciptype(nv,ifunit,jf,iens,
     +      ptype_mn,ptype_pr)

            use grib_mod
            include 'parm.inc'

  
c    for derived variables
        Character*4 dvname(maxvar)
        Integer dk5(maxvar), dk6(maxvar),dk4(maxvar)
        Character*1 dMsignal(maxvar), dPsignal(maxvar)
        Integer dMlvl(maxvar), dMeanLevel(maxvar,maxmlvl)
        Integer dPlvl(maxvar), dProbLevel(maxvar,maxplvl)
        Character*1 dop(maxvar)
        Integer dTlvl(maxvar)
        Real    dThrs(maxvar,maxtlvl)
        Integer MPairLevel(maxvar,maxmlvl,2)
        Integer PPairLevel(maxvar,maxplvl,2)
                                                                                                                                                                
        common /dtbl/nderiv,
     +              dvname,dk4,dk5,dk6,dMlvl,dPlvl,dTlvl,
     +              dMeanLevel,dProbLevel,dThrs,
     +              dMsignal,dPsignal,MPairLevel,PPairLevel,dop

        INTEGER, intent(IN) :: nv, jf, iens

        REAL,dimension(jf,4),intent(INOUT) :: ptype_mn
        REAL,dimension(jf,4),intent(INOUT) :: ptype_pr

        INTEGER miss(iens)
        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld 
        integer  jpd1(4),jpd2(4),jpd10(4),jpd12(4)
        real prcptype(jf,iens,4)                 !store 4 precip types data

        data (jpd1(i),i=1,4)
     +  /1,1,1,1/
        data (jpd2(i),i=1,4)
     +  /192,193,194,195/                     !rain, frzr, icep, snow
        data (jpd10(i),i=1,4)
     +  /1,1,1,1/
        data (jpd12(i),i=1,4)
     +  /0,0,0,0/

        !write(*,*) ' In preciptype ......'
        !write(*,*) nv,ifunit,jf,iens

        miss=0
        mbrs=iens
        !First get necessary data from raw ensemble members 
        jpd27=-9999
        jpdtn=0
        DO 101 irun=1,iens
          do i=1,4
           call readGB2(ifunit(irun),jpdtn,jpd1(i),jpd2(i),
     +       jpd10(i),jpd12(i),jp27,gfld,iret)
           if(iret.eq.0) then
             prcptype(:,irun,i)=gfld%fld
           else
             miss(irun)=1
             mbrs=mbrs - 1
             goto 101
           end if
          end do
101     CONTINUE

c         write(*,*) 'Get data done '

           ptype_mn = 0.
           ptype_pr = 0.

              
          do 200 lv=1,dMlvl(nv)                                  !dMlvl(nv) always =1 for precp type
            do 300 igrid = 1,jf
         
                 crain = 0.
                 cfrzr = 0.
                 csnow = 0.
                 cslet = 0.

                 do irun = 1,iens 
                  if(miss(irun).eq.0) then
                   crain = crain + prcptype(igrid,irun,1)
                   cfrzr = cfrzr + prcptype(igrid,irun,2)
                   cslet = cslet + prcptype(igrid,irun,3) 
                   csnow = csnow + prcptype(igrid,irun,4)
                  end if
                 end do  

cc  following is part is the code copy from Geoff Manikin doninant precip type decision
cc  importance priority order:  freezing_rain(1) > snow(2) > sleet(3) > rain(4)            
          if(crain.ge.1.0.or.cfrzr.ge.1.0.or.csnow.ge.1.0.or.
     +         cslet.ge.1.0) then
            if(csnow.ge.cslet) then
              if(csnow.ge.cfrzr) then
                  if(csnow.ge.crain) then
                   ptype_mn(igrid,4)=1.      !snow
                   goto 800
                  else
                   ptype_mn(igrid,1)=1.      !rain
                   goto 800
                  end if
                else if(cfrzr.ge.crain) then
                   ptype_mn(igrid,2)=1.      !freezing rain
                   goto 800
                else
                   ptype_mn(igrid,1)=1.      !rain
                   goto 800
                end if
             else if(cslet.gt.cfrzr) then
                if(cslet.ge.crain) then
                   ptype_mn(igrid,3)=1.      !sleet
                   goto 800
                else
                   ptype_mn(igrid,1)=1.      !rain
                   goto 800
                end if
             else if(cfrzr.ge.crain) then
                   ptype_mn(igrid,2)=1.      !freezing rain
                   goto 800
             else
                   ptype_mn(igrid,1)=1.      !rain
                   goto 800
               end if

800          continue    
            end if

c             if(igrid.eq.150000) then
c              write(*,*)'IN preciptype 150000 ptype_mn'
c              write(*,*)cfrzr,csnow,cslet,crain,
c     +        (ptype_mn(igrid,n),n=1,4)
c             end if

300        continue 
200      continue 

        do 400 lv=1,dPlvl(nv)     !dPlvl(nv) always is 1 for precp type
          do 500 lt =1,dTlvl(nv)
            do 600 igrid = 1, jf
   
                 crain = 0.
                 cfrzr = 0.
                 csnow = 0.
                 cslet = 0.

                if(dThrs(nv,lt).eq.1.0) then

                 do irun=1,iens
                   if(miss(irun).eq.0) then
                   crain = crain + prcptype(igrid,irun,1)
                   cfrzr = cfrzr + prcptype(igrid,irun,2)
                   cslet = cslet + prcptype(igrid,irun,3)
                   csnow = csnow + prcptype(igrid,irun,4) 
                  end if
                 end do
             
                  if (mbrs.gt.0) then
                     ptype_pr(igrid,1)=crain*100./mbrs
                     ptype_pr(igrid,2)=cfrzr*100./mbrs        
                     ptype_pr(igrid,3)=cslet*100./mbrs
                     ptype_pr(igrid,4)=csnow*100./mbrs
                  else
                     ptype_pr(igrid,1)=0.
                     ptype_pr(igrid,2)=0.
                     ptype_pr(igrid,3)=0.
                     ptype_pr(igrid,4)=0.
                  end if   
                end if


c             if(igrid.eq.150000) then
c              write(*,*)'IN preciptype 150000 ptype_pr'
c              write(*,*) (ptype_pr(igrid,kp),kp=1,4)
c             end if

600      continue
500     continue
400    continue

           return
           end
