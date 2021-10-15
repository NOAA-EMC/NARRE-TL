cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine flight_res: compute flight_restriction condition prob
c     
c     Author: Binbin Zhou, Oct 1, 2007
c     Modification: 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine flight_res (nv,ifunit,jf,iens,Lp,Lt,
     +              derv_pr,wgt)

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
        REAL,dimension(jf,Lp,Lt),intent(INOUT) ::  derv_pr

        
        real  count, aprob, flt_cnd(iens),wgt(30),visb,ceil
        integer ID_FLT

        integer miss(iens)

        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld

        real,dimension(jf,iens) :: ceiling,visbil

        jpdtn=0
        jp27=-9999

        write(*,*) 'In flight_res .....'
        write(*,*) 'nv,ifunit,jf,iens,Lp,Lt,jpd10',
     +              nv,ifunit,jf,iens,Lp,Lt,jpd10

        miss=0
    
        do 200 k=1,iens

          jpd12=dProbLevel(nv,1)

           call readGB2(ifunit(k),jpdtn,3,5,215,0,jp27,gfld,iret)  !Sfc height
            if(iret.eq.0) then
             ceiling(:,k)=gfld%fld
            else
             miss(k)=1
             !write(*,*) 'Sfc height  missing in file',ifunit(k) 
             goto 200
            end if

           call readGB2(ifunit(k),jpdtn,19,0,1,0,jp27,gfld,iret) !Sfc visb
            if(iret.eq.0) then
             visbil(:,k)=gfld%fld
            else
             miss(k)=1
             !write(*,*) ' Sfc visb  missing in file',ifunit(k)
             goto 200
            end if

 200   continue  

        do 600 igrid=1,jf
           do k=1,iens
            if(miss(k).eq.0) then
              ceil=ceiling(igrid,k)
              visb=visbil(igrid,k)

              call flight_cond(ceil,visb,fltc)

              !if(igrid.ge.10000.and.igrid.le.10100) then
              !  write(*,*) igrid,k,tcld,cldb,sfch,visb,fltc
              !end if 

               flt_cnd(k)=fltc
             else
               flt_cnd(k)=0.0
             end if
           end do
 

           !do k=1,iens                       !this is to exclude spc's wrong vis data
           ! visb=visbil(igrid,k)
           ! if(visb.lt.1.0) miss(k)=1
           !end do
        
           do lv=1,dPlvl(nv)
            do lt = 1, dTlvl(nv)
 
             thr1 = dThrs(nv,lt)
             thr2 = 0.
             call getprob(flt_cnd,iens,thr1,thr2,dop(nv),aprob,
     +         miss,wgt)
              
             derv_pr(igrid,lv,lt)=aprob
           end do
          end do

600      continue

        return
        end


C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .
C SUBPROGRAM:    FLIGHT_CONDD   COMPUTES flight condition
C   PRGRMMR: Binbin Zhou      /NCEP/EMC  DATE: 2005-08-18
C
C ABSTRACT:
C    This program computes the flight condition restriction
C    which is defined as follow (NOAA/NWS/Instruction for TAF, 2004):
C 
C                Ceiling(feet)             Visibility(miles)   FLTCND
C      LIFR        < 200           and/or      < 1               1
C      IFR      >= 500 to <  1000  and/or     >=1 to <  3        2
C      MVFR     >=1000 to <= 3000  and/or     >=3 to <= 5        3
C      VFR         > 3000                       > 5              4
C


        subroutine flight_cond(ceil,visb,fltc)

         real visb,fltc,ceil

          CEIL=CEIL*3.2808      !m->feet
          VISB = VISB / 1609.0  !m-> mile

         !compute flight condition

         IF(CEIL.LT.500.0 .OR. VISB.LT.1.0 ) THEN
             FLTC = 1.0

          ELSE IF( (CEIL.GE.500.AND.CEIL.LT.1000.0) .OR.
     +              (VISB.GE.1.0.AND.VISB.LT.3.0) ) THEN
             FLTC = 2.0

          ELSE IF( (CEIL.GE.1000.AND.CEIL.LE.3000.0) .OR.
     +              (VISB.GE.3.0.AND.VISB.LE.5.0) ) THEN
             FLTC = 3.0

          ELSE IF( CEIL.GT.3000.0  .OR. VISB.GT.5.0) THEN
             FLTC = 4.0

          END IF

         return
         end

