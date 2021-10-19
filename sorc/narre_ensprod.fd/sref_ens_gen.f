C$$$me  MAIN PROGRAM DOCUMENTATION BLOCK $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C MAIN PROGRAM: SREF_ENS_GEN
C   PRGMMR: DU               ORG: NP21        DATE: 2002-01-02
C
C ABSTRACT: READS IN SREF GRIB AND OUTPUTS ENSEMBLE GRIB 
C           PRODUCTS. 
C
C PROGRAM HISTORY LOG:
c
c Old Version:
c 2001-07-10    Jun Du, initial program
c
c Generic Version: 
c
c 2005-08-01 Author: Binbin Zhou 
c
c       Re-organize the code structure and re-write the program in following features:
c
c        (1) The generic version dosn't fix GRIB# like 212 as old version. It aan process
c            any GRIB#, and the domain (IM, JM) for the defined GRIB# is also flexible. 
c            As long as GRIB# is specified in 'filename', this program can automatically 
c            retrievs the domain (IM, JM, JF=IM*JM), then all arrays associated with this 
c            domain are dynamically allocated (on the fly)
c        (2) Ensemble members (IENS) are not fixed, but depends on the input file "filename'
c            in which, the all of soft linked file name are listed in order. The order and
c            number of files are changeable, the program can automatically allocate them
c        (3) The variables are changeable by using the 'variable.tbl' file, in which, two kinds of
c            variables are defined: (i) direct variables: read from GRIB files. Direct  
c            variables will also be computed for mean/spread an probability depending 
c            on the settings in each variable redcord (ii) derived variables: they are 
c            not read from GRIB files but need to do some derivation from direct variables.
c        (4) Derived variables need some algorithm to derive, Currently, only following derived 
c            variables are defined in this program:
c            (A) Thickness of two (pressure) levels: mean/spread/probability         
c            (B) 3hr, 6hr, 12hr, 24hr Accumulated precipitation: mean/spread/probability 
c            (C) 3hr, 6hr, 12hr, 24hr Accumulated snowfall: mean/spread/probability
c            (D) Precipitation type: probability
c            (E) Dominant precipitation: Mean
c            If want to process other derived product, user must provide their algorithms 
c            and then the computation code should be put into this program       
c        (5) User can request any kind of direct variables, either mean/spread or probability
c            by editing the 'variable.tbl' file (read README file for its gramma)
c 
c change log:
c 2006-01-12:Binbin Zhou
c             Add Geoff Manikin's algorithm to determine dominant precip type 
c 2006-03-01:Binbin Zhou
c             Modify wind spread computation by considering the wind directions, suggested 
c             by Jun Du            
c 2006-03-30: Jun Du: convert the unit from K to C for lifted index 
c
c 2006-04-30: Binbin Z. Add DTRA variance variables 
c
c 2006-05-15: Binbin Zhou
c             Modify to accept Beijing 2008 Olympic Game Grid#137
c
c 2006-07-21: Binbin Z. Add fog probability
c
c 2007-01-09: Binbin Z. Add freezing rain for GFS ensemble
c
c 2007-01-15: Binbin Z. Add Absolute vorticity for global ensemble since GFS has no absolute vorticity output
c 
c 2007-04-05: Binbin Z. Add High resolution WRF testing grid#137 
c
c 2007-10-09: Binbin Z. Add flight-restriction probability computation
c 
c 2009-04-07: Binbin Z. Adopted for VSREF (by adding member weights and RUC/VSREF's grid#130 ) 
c
c 2009-09-09: Binbin Z. Add convection adapted from Steve W. of GSD
C
c 2010-06-20: Binbin Z. Add precip accumulation for VRSEF version
c
c 2010-8-18:  Binbin Z. Add new fog algorithm (based on Zhou and Ferrier 2008)
c 2011-3-18:  Binbin Z. adapt for NARRE-Time Lag  (based on Zhou and Ferrier 2008)
c
c 2011-04-15: Binbin Z. Add cpat lightning products (based on David Bright's SPC program)
c 2011-04-18: Binbin Z. Add fire weather probability computation
c 2011-04-18: Binbin Z. Add LLWS code
c 2011-06-10: Binbin Z. Add cptp_dryt based on David B.'s code (but dryt_hrly_rgn3 seems not good)
c 2011-06-20: Binbin Z. Add severe thunder storm potential based on David B's code
c 2012-01-20: Binbin Z. Add missing array to deal with after copygb of hiresw WRF-east/west over CONUS 
c                       Modified to work on NSSE 
c 2012-02-21: Binbin Z. Modify accumulated precip computation method: read them directly from prcip.* files
c                       so they become direct vriables in the table 
c 
c 2013-03-28: Binbin Z. Modify to work on grib2 and structure change to reduce memory by read one variable
c                       then process it and pack its ensmeble products into output files
c
c 2013-12-21: Binbin Z. Modify  missing array from missing(jf,iens) to missing(maxvar,iens) after 
c                       hiresWRFs are all in same CONUS grid, and deal with field missing in GRIB2 file 
c                       and getGB2 issue (not stop if the field is not in the GRIB2 file
c
c 2014-03-15: Binbin Z. Add Storm-related fileds requested by Jacob Caley of EMC
c
c 2014-04-21: Binbin Z. Add Flash Flood and Intense Rain (FFaIR) summer experiment fields
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C$$$
C

       use grib_mod
       include 'parm.inc'

      type(gribfield) :: gfld, gfld_temp, gfld_nam,gfld_rap,gfld_vis,
     +                   gfld_ceil
C  raw data
      real,allocatable,dimension(:,:,:)   :: rawdata_mn, rawdata_pr   !jf,iens, maxmlvl
      real,allocatable,dimension(:,:,:)   :: precip                   !jf,iens, number of fcst output files
      real,allocatable,dimension(:,:,:)   :: mrk_ice                  !jf,iens, number of fcst output files
      real,allocatable,dimension(:,:,:)   :: mrk_frz                  !jf,iens, number of fcst output files
      real,allocatable,dimension(:,:,:)   :: mrk_snow                 !jf,iens, number of fcst output files

C mean
      real,allocatable,dimension(:,:) :: vrbl_mn                    !jf, maxmlvl
      real,allocatable,dimension(:,:) :: derv_mn                    !jf, maxmlvl

C spread
      real,allocatable,dimension(:,:) :: vrbl_sp                    !jf, maxmlvl
      real,allocatable,dimension(:,:) :: derv_sp                    !jf, maxmlvl

      real,allocatable,dimension(:) :: ppt3_sp,ppt6_sp,ppt12_sp,
     *                                 ppt24_sp,s12_sp

C probability
      real,allocatable,dimension(:,:,:) :: vrbl_pr                   !jf, maxplvl, maxtlvl
      real,allocatable,dimension(:,:,:) :: derv_pr                   !jf, maxplvl, maxtlvl

C Mix,min, 10,25,50 and 90% mean products
      real,allocatable,dimension(:,:,:) :: mxp8                     !jf,maxmlvl,7


C others 
       character*19, allocatable, dimension(:) :: fhead
       real,allocatable,dimension(:) ::           p03mp01              !jf
       real,allocatable,dimension(:,:)  :: ptype_mn,ptype_pr           !jf,4 

       real,allocatable,dimension(:,:,:) :: derv_dtra                  !jf,maxmlvl,8 for DTRA requests
       real,allocatable, dimension(:)    :: Hsfc                       !surface height for DTRA  
 
       integer,allocatable,dimension(:,:)   :: missing                 ! to deal with missing data 
       integer,allocatable,dimension(:)     :: miss
       real,allocatable,dimension(:) ::           apoint               !iens

C for get grib size jf=im*jm
  
       integer iyr,imon,idy,ihr
       character*50 gdss(400)
       integer IENS, gribid, kgdss(200), lengds,im,jm,km,jf

       integer  leadtime, interval, loutput

cc%%%%%%%  8. To be added if necessary ...........................
C original                                                                                                                                           
       dimension kgds(25)
       character*20 mnout,pmin,pmax,pmod,pp10,pp25,pp50,pp75,pp90
       character*20 spout
       character*20 prout
       character*40 files(50),prcps(50)
     


C for variable table:
        Integer numvar, nderiv
        Character*4 vname(maxvar)
        Integer k5(maxvar), k6(maxvar), k4(maxvar)     !k4-category#, k5-paremter#, k6-surface id
        Character*1 Msignal(maxvar), Psignal(maxvar)
        Integer Mlvl(maxvar), MeanLevel(maxvar,maxmlvl),Lm
        Integer Plvl(maxvar), ProbLevel(maxvar,maxplvl),Lp
        Integer Tlvl(maxvar),Lth
        Character*1 op(maxvar)
        Real    Thrs(maxvar,maxtlvl)
        Character*5 eps       
  
c    for derived variables
        Character*4 dvname(maxvar)
        Integer dk5(maxvar), dk6(maxvar), dk4(maxvar)
        Character*1 dMsignal(maxvar), dPsignal(maxvar)
        Integer dMlvl(maxvar), dMeanLevel(maxvar,maxmlvl)
        Integer dPlvl(maxvar), dProbLevel(maxvar,maxplvl)
        Character*1 dop(maxvar)
        Integer dTlvl(maxvar)
        Real    dThrs(maxvar,maxtlvl)
        Integer MPairLevel(maxvar,maxmlvl,2)
        Integer PPairLevel(maxvar,maxplvl,2)

c   for max,min,10,25,50,90% mean products
        Character*4 qvname(maxvar)
        Integer qk5(maxvar), qk6(maxvar), qk4(maxvar)
        Character*1 qMsignal(maxvar)
        Integer qMlvl(maxvar), qMeanLevel(maxvar,maxmlvl)
        
        real  weight(30)                               
        character*20 filenames
        character*2 cfhr                                       
        integer est                                 !east time for convection code
        integer ifunit(50),ipunit(50)               !member, prcp member files units         
        character*2 cycle(24)
        character*7 mbrname(50)

        integer jptyp2(4)                           !for precip type jpd2
        integer iqout(8) 

        integer ierr_open_punit(30),ierr_open_funit(30),
     +          ierr_open(30) 

        common /tbl/numvar,
     +              vname,k4,k5,k6,Mlvl,Plvl,Tlvl,
     +              MeanLevel,ProbLevel,Thrs,
     +              Msignal,Psignal,op
                                                                                                                                                                
        common /dtbl/nderiv,
     +              dvname,dk4,dk5,dk6,dMlvl,dPlvl,dTlvl,
     +              dMeanLevel,dProbLevel,dThrs,
     +              dMsignal,dPsignal,MPairLevel,PPairLevel,dop
                           
     
        common /qtbl/nmxp,
     +              qvname,qk4,qk5,qk6,qMlvl,
     +              qMeanLevel,qMsignal
                                                                                                                                
        data (iqout(i),i=1,8)
     +  /301,302,303,304,305,306,307,308/

        data (cycle(i),i=1,24)
     +  /'00','01','02','03','04','05','06','07','08',
     +   '09','10','11','12','13','14','15','16','17',
     +   '18','19','20','21','22','23'/

        data (jptyp2(i),i=1,4)
     +  /192,193,194,195/                     !rain, frzr, icep, snow

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c (I)   Read table file to get wanted ensemble variable information
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        numvar=0
        nderiv=0
        nmxp=0
        ierr_open_punit=0
        ierr_open_funit=0
        ierr_open=0

        nunit=101
        open (nunit, file='variable.tbl', status='old')
        
        write(*,*) 'reading variable.tbl ...'                           
        call readtbl(nunit)
        write(*,*) 'read variable.tbl done'

        close(nunit)

c        write(*,*) 'Check rirect variable reading:'
c        do i = 1, numvar
c          write(*,*) vname(i),k4(i),k5(i), k6(i),Msignal(i), Mlvl(i),
c     +                 (MeanLevel(i,j),j=1,Mlvl(i)),
c     +    Psignal(i), Plvl(i), (ProbLevel(i,j),j=1,Plvl(i)),
c     +    op(i), Tlvl(i),   (Thrs(i,j),j=1,Tlvl(i))
c151      format(a4,1x,3i4,a2,i2,<Mlvl(i)>(1x,i4),a2,i2,<Plvl(i)>(1x,i4),
c     +   a2, i2, <Tlvl(i)>(1x,F7.1))   
c        end do

c        write(*,*) 'Check derived variable reading:'
c        do i = 1, nderiv
c         if(dk6(i).eq.101) then
c           write(*,*) dvname(i),dk4(i),dk5(i), dk6(i),dMsignal(i),
c     +               dMlvl(i),(dMeanLevel(i,j),j=1,dMlvl(i)),
c     +              (MPairLevel(i,j,1),j=1,dMlvl(i)),
c     +              (MPairLevel(i,j,2),j=1,dMlvl(i)),
c     +          dPsignal(i),dPlvl(i),(dProbLevel(i,j),j=1,dPlvl(i)),
c     +              (PPairLevel(i,j,1),j=1,dPlvl(i)),
c     +              (PPairLevel(i,j,2),j=1,dPlvl(i)),
c     +              dop(i), dTlvl(i),   (dThrs(i,j),j=1,dTlvl(i))
c         else
c           write(*,*) dvname(i),dk4(i),dk5(i), dk6(i),dMsignal(i),
c     +              dMlvl(i),(dMeanLevel(i,j),j=1,dMlvl(i)),
c     +          dPsignal(i),dPlvl(i),(dProbLevel(i,j),j=1,dPlvl(i)),
c     +             dop(i), dTlvl(i),   (dThrs(i,j),j=1,dTlvl(i))
c         end if
c        end do
c
c        write(*,*) 'Check max,min,10,25,50,75,90% mean product reading:'
c        do i = 1, nmxp
c         write(*,*)qvname(i),qk4(i),qk5(i),qk6(i),qMsignal(i),
c     +        qMlvl(i),(qMeanLevel(i,j),j=1,qMlvl(i))
c        end do


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c (II): Read filename file to get # of ensemble available, grib#, vertical pressure level,
c       last forecast time, forecast outout interval and grib file heads, 
c       so that arrays can be dinamically allocatable
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       filenames="filename" 
       write(*,*) 'read ',filenames
       open(1, file=filenames, status='old')
       read(1,*) iyr,imon,idy,ihr,ifhr,gribid,KM,leadtime,interval, !ihr-cycle hour , ifhr-forecast hour
     +             loutput
       write(*,*)iyr,imon,idy,ihr,ifhr,gribid,KM,leadtime,interval,
     +             loutput

       !read(1,*) (cfhr(k),k=1,loutput) 

       IENS=0
       do irun=1,50
         read(1,301,END=302) weight(irun),files(irun),mbrname(irun)
         write(*,301) weight(irun), files(irun),mbrname(irun)
         IENS=IENS+1
         prcps(irun)='prcip'//files(irun)(5:17)              !get prcip file names
         write(*,*) weight(irun), files(irun),prcps(irun)
       end do

       close (1)

  301  format(f7.2, 1x, a17,4x,a7)
  302  continue
      write(*,*) 'IENS=',IENS, (weight(irun),irun=1,iens)

      DO 3000 itime=ifhr,ifhr          !Here ifhr is forecast hour, e.g. .f12, '12' is ifhr

         write(cfhr,'(i2.2)') itime

        print*,'Initializing variables'

       est = ihr + itime - 6                                          !Eastern time For convection computation
       if (est.lt.0) est= est + 24 


CCCC Binbin Zhou Note:

       if(gribid.eq.255) then   !For HRRR grid
         im=1799
         jm=1059
         jf=im*jm
         dx=3000.0
         dy=3000.0
       else
         call makgds(gribid, kgdss, gdss, lengds, ier)
         im=kgdss(2)
         jm=kgdss(3)
         jf=kgdss(2)*kgdss(3)
         dx=1.*kgdss(8)
         dy=1.*kgdss(9)
       end if
 
       write(*,*) 'im, jm, jf,dx,dy =', im, jm, jf,dx,dy
       write(*,*) 'leadtime, interval, loutput=', 
     +             leadtime, interval, loutput

       if (.NOT.allocated(fhead)) then
        allocate(fhead(iens))   
       end if
       if (.NOT.allocated(apoint)) then
        allocate(apoint(iens)) 
       end if
       if (.NOT.allocated(missing)) then
        allocate(missing(maxvar,iens))
       end if
       if (.NOT.allocated(miss)) then
        allocate(miss(iens))
       end if

       if (.NOT.allocated(p03mp01)) then
        allocate(p03mp01(jf))
       end if



       eps=files(1)(1:4)
       missing=0 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c (III): Open all ensemble data files 
c In these lines put here, we need the date and time of the forecast,
c the name of the file.  We need to open the file.  Then we can point to
c the correct location of the forecast in the forecast file.  If this
c is an ensemble, then we also need another loop for the forecast name
c (cntl, n1, p1, etc.).  That loop should be inside the date/time loop.


       do 500 irun=1,iens

        ifunit(irun)=10+irun
        call baopenr(ifunit(irun),files(irun),ier1)
        write(*,*)'open#',irun,ifunit(irun),files(irun),'err=',ier1 
        if (ier1.ne.0) ierr_open_funit(irun)=ier1

        ipunit(irun)=50+irun
        call baopenr(ipunit(irun),prcps(irun),ier2)
        write(*,*)'open#',irun,ipunit(irun),prcps(irun),'err=',ier2
        if (ier2.ne.0) ierr_open_punit(irun)=ier2

500     continue 

        write(*,*) 'numvar=',numvar

        imean=201
        isprd=202
        iprob=203

      mnout=trim(eps)//'.mean.t'//cycle(ihr+1)//'z'//'.f'//cfhr
      spout=trim(eps)//'.sprd.t'//cycle(ihr+1)//'z'//'.f'//cfhr
      prout=trim(eps)//'.prob.t'//cycle(ihr+1)//'z'//'.f'//cfhr

      call baopen (imean, mnout, iret)
      if (iret.ne.0) write(*,*) 'open ', mnout, 'err=', iret
      call baopen (isprd, spout, iret)
      if (iret.ne.0) write(*,*) 'open ', spout, 'err=', iret
      call baopen (iprob, prout, iret)
      if (iret.ne.0) write(*,*) 'open ', prout, 'err=', iret


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c (IV) Main tasks:
c  Loop1: process direct variables
c  Loop2: process derived variables
c  Loop3: process mxp variables (min,max,mod, 10, 250 50 75 90%) 

c  Loop 1: for direct variables
c   1-0: Allocate nesessary arrays
c   1-1: In this loop, one field (nv) data are read from all members (1001 loop),
c   1-2: Compute its mean, spread, and prob,
c   1-3: Pack them into output mean, spread, prob files   
c   then, loop for next field
c   1-4: deallocate allocated arrays 
c
       DO 2001 nv = 1, numvar

c       write(*,*) 'nv Mlvl(nv),Plvl(nv),Tlvl(nv)=',
c     +    nv, Mlvl(nv),Plvl(nv),Tlvl(nv)

c  Loop  1-0: Allocate nesessary arrays
         if (.NOT.allocated(rawdata_mn)) then
           allocate (rawdata_mn(jf,iens,Mlvl(nv)))
         end if
         if (.NOT.allocated(rawdata_pr)) then
           allocate (rawdata_pr(jf,iens,plvl(nv)))
         end if
        if (.NOT.allocated(vrbl_mn)) then
           Lm=max(1,Mlvl(nv))                !Keep at least one vrbl_mn to pass it into packGB2 
           allocate (vrbl_mn(jf,Lm))         !in case Mlvl(nv)=0
         end if
        if (.NOT.allocated(vrbl_sp)) then
           Lm=max(1,Mlvl(nv))
           allocate (vrbl_sp(jf,Lm))
         end if
        if (.NOT.allocated(vrbl_pr)) then
           Lp=max(1,Plvl(nv))
           Lth=max(1,Tlvl(nv))
           allocate (vrbl_pr(jf,Lp,Lth))
         end if


         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Loop 1-1: Read direct variable's GRIB2 data from all members 



        DO 1001 irun=1,iens         ! members

         !apcp1h,apcp3hr,apcp6hr,apcp12,apcp24 have been put into 
         !every members as direct variable

          if ((k4(nv).eq.2.and.k5(nv).eq.222).or.
     +        (k4(nv).eq.2.and.k5(nv).eq.223).or.
     +        (k4(nv).eq.16.and.k5(nv).eq.198).or.
     +        (k4(nv).eq.16.and.k5(nv).eq.195).or.
     +        (k4(nv).eq.16.and.k5(nv).eq.196).or.
     +        (k4(nv).eq.7.and.k5(nv).eq.199).or.
     +        (k4(nv).eq.2.and.k5(nv).eq.220).or.
     +        (k4(nv).eq.2.and.k5(nv).eq.221) ) then

              if(mbrname(irun).eq.'hrrrgsd') then ! These 8 fields are not in GSD's HRRR
                 rawdata_pr(:,irun,:)=-9999.0
                 rawdata_mn(:,irun,:)=-9999.0
                 goto 1001
              end if

           end if  

          if((k4(nv).eq.1.and.k5(nv).eq.8).or.          !Check if var is precip/snow to set file unit 
     +      (k4(nv).eq.1.and.k5(nv).eq.15)) then

c         write(*,*) 'get APCP GRIB2 data for member ', irun

           jpdtn=8
           jpd1=k4(nv)
           jpd2=k5(nv)
           jpd10=k6(nv)
           !jpd12 is determined by a specific level MeanLevel(nv,lv) to !ProbLevel(nv,lv) later on

           if (vname(nv).eq.'AP1h'.or.vname(nv).eq.'SN1h') then        !AP1h,Ap3h,Ap6h, AP12,Ap24 should be hardcopy in the variable tbl
             jpd27=1
           else if (vname(nv).eq.'AP3h'.or.vname(nv).eq.'SN3h') then
             if(ifhr.lt.3) goto 222
             jpd27=3
           else if (vname(nv).eq.'AP6h'.or.vname(nv).eq.'SN6h') then
             if(ifhr.lt.6) goto 222
             jpd27=6
           else if (vname(nv).eq.'AP12'.or.vname(nv).eq.'SN12' ) then
             if(ifhr.lt.12) goto 222
             jpd27=12
           else if (vname(nv).eq.'AP24'.or.vname(nv).eq.'SN24' ) then
             if(ifhr.lt.24) goto 222
             jpd27=24
           else
             write(*,*) 'Using wrong APCP name!'
           end if

           igrb2=ipunit(irun)
           ierr_open(irun)=ierr_open_punit(irun)
          else  

           jpdtn=0

          if ((k4(nv).eq.2.and.k5(nv).eq.222).or.
     +        (k4(nv).eq.2.and.k5(nv).eq.223).or.
     +        (k4(nv).eq.16.and.k5(nv).eq.198).or.
     +        (k4(nv).eq.7.and.k5(nv).eq.199).or.
     +        (k4(nv).eq.2.and.k5(nv).eq.220).or.
     +        (k4(nv).eq.2.and.k5(nv).eq.221) ) then
              jpdtn=8
           end if

           jpd1=k4(nv)
           jpd2=k5(nv)
           !jpd12 is determined by a specific level MeanLevel(nv,lv) to !ProbLevel(nv,lv) later on
           jpd10=k6(nv)
           jpd27=-9999

           igrb2=ifunit(irun)
           ierr_open(irun)=ierr_open_funit(irun)

          end if

         do 501 lv=1,Mlvl(nv)

          jpd12=MeanLevel(nv,lv)

          if(mbrname(irun).eq.'hrrrgsd'.and.k4(nv).eq.16.
     +                            and.k5(nv).eq.197) then !HRRR use different echo-top ID from standard table
            jpd1=3
            jpd2=193
            jpd10=3
            jpd12=-9999
          end if  

          if (ierr_open(irun).eq.0) then
            call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +            gfld, jret)
            if (jret.ne.0) goto 501
            rawdata_mn(:,irun,lv)=gfld%fld
            !if (nv.eq.1.and.irun.eq.8) gfld_nam=gfld  !NAM's gfld, order of files in filename must be: first 7 are rap, last 4 are nam
            if (nv.eq.1.and.irun.eq.5) gfld_rap=gfld  !RAP's gfld, order of files in filename must be: first 7 are rap, last 4 are nam
            !if (nv.eq.1.and.irun.eq.1) gfld_rap=gfld  !RAP's gfld
            if (nv.eq.1.and.irun.eq.1) gfld_nam=gfld  !NAM's gfld
            if (k4(nv).eq.19.and.k5(nv).eq.0.and.k6(nv).eq.1   !visibility mean uses RAP's gfld to packGB2_mean
     +          .and.irun.eq.8)  gfld_vis=gfld
            if (k4(nv).eq.3.and.k5(nv).eq.5                    !ceiling mean uses RAP's gfld to packGB2_mean
     +          .and.irun.eq.8)  gfld_ceil=gfld

          end if
501      continue

         do 502 lv=1,Plvl(nv)

          jpd12=probLevel(nv,lv)

           if(mbrname(irun).eq.'hrrrgsd'.and.k4(nv).eq.16.
     +                            and.k5(nv).eq.197) then !HRRR use different echo-top ID from standard table
            jpd1=3
            jpd2=193
            jpd10=3
            jpd12=-9999
           end if

           if (ierr_open(irun).eq.0) then
            call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +          gfld,kret)
            if (kret.eq.0) then
             rawdata_pr(:,irun,lv)=gfld%fld
            else 
              if (jpd1.eq.16.and.jpd2.eq.196        !RAP REFC jpd10=10 instead of 200
     +           .and.jpd10.eq.200) then
                write(*,*) 'Try jpd10=10 to read REFC RAP',irun
                call readGB2(igrb2,jpdtn,jpd1,jpd2,10,jpd12,jpd27,
     +          gfld,kret)
               if (kret.eq.0) rawdata_pr(:,irun,lv)=gfld%fld
              end if        
            end if
          end if 

502       continue       
 
         !Just correc some variables 
         if(jret.ne.0.or.kret.ne.0) then
           missing(nv,irun)=1
           goto 1001
         end if 

         do igrid=1,jf
          do lv = 1, Mlvl(nv)

            if(k4(nv).eq.7.and.k5(nv).eq.192) then
              rawdata_mn(igrid,irun,lv) = rawdata_mn(igrid,irun,lv)
     +                                                       - 273.15  !Jun Du: change unit from K to C for lifted index
            end if

            if(k4(nv).eq.3.and.k5(nv).eq.5.and.(k6(nv).eq.2.or.
     +        k6(nv).eq.3))                                    then    !RUC cloud base/top with negative values

              if (rawdata_mn(igrid,irun,lv).le.0.0) then               !means no cloud
                 rawdata_mn(igrid,irun,lv) = 20000.0
              end if
            end if

           end do
 
           do lv = 1, Plvl(nv)
             if(k4(nv).eq.7.and.k5(nv).eq.192) then
               rawdata_pr(igrid,irun,lv) = rawdata_pr(igrid,irun,lv)
     +                                                       - 273.15  !Jun Du: change unit from K to C for lifted index
             end if

            if(k4(nv).eq.3.and.k5(nv).eq.5.and.(k6(nv).eq.2.or.
     +        k6(nv).eq.3))                                    then     !RUC cloud base/top with negative values

              if (rawdata_pr(igrid,irun,lv).le.0.0) then               !means no cloud
                rawdata_pr(igrid,irun,lv) = 20000.0
              end if
            end if
           end do

         end do

1001     continue  !end of iens loop for getting GRIB2 data for this direct variable
       
         write(*,*) 'Get direct variable data done for', nv

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Loop 1-2:   Compute mean/spread/prob for this direct variable

         vrbl_mn=0.
         vrbl_sp=0.
         vrbl_pr=0.

        IF(trim(Msignal(nv)).eq.'M') THEN
         
         do lv = 1, Mlvl(nv)
          do igrid=1,jf

           apoint = rawdata_mn(igrid,:,lv)
           miss = missing(nv,:)

          if ((k4(nv).eq.2.and.k5(nv).eq.222).or.   !These field not in hrrrgsd file
     +        (k4(nv).eq.2.and.k5(nv).eq.223).or.
     +        (k4(nv).eq.16.and.k5(nv).eq.198).or.
     +        (k4(nv).eq.16.and.k5(nv).eq.195).or.
     +        (k4(nv).eq.16.and.k5(nv).eq.196).or.
     +        (k4(nv).eq.7.and.k5(nv).eq.199).or.
     +        (k4(nv).eq.2.and.k5(nv).eq.220).or.
     +        (k4(nv).eq.2.and.k5(nv).eq.221) ) then
             do ie=1,iens
               if(mbrname(ie).eq.'hrrrgsd') then
                miss(ie)=1
               end if
             end do
          end if


           if(k4(nv).eq.19.and.k5(nv).eq.0) then  ! Visibility: conditional mean 
             do ie =1,iens
               if(apoint(ie).lt.1.0) miss(ie)=1   !exclude spc runs's error vis 
             end do

              
             call get_cond_mean (apoint,iens,24056.0,amean,aspread,
     +                              miss, weight) 

           else if(k4(nv).eq.3.and.k5(nv).eq.5.and.k6(nv).eq.215)then !  Ceiling: conditional mean
             call get_cond_mean (apoint,iens,20000.0,amean,aspread,
     +                               miss,weight)

           else if(k4(nv).eq.6.and.k5(nv).eq.6.and.k6(nv).eq.1) then  !  Sfc CLOUD (fog): liquid water content
                apoint = apoint * 1000.0                              !  kg/kg -> g/kg
                call get_cond_mean_lwc(apoint,iens,0.0,amean,aspread,
     +                               miss,weight)

           else if(k4(nv).eq.3.and.k5(nv).eq.5.and.k6(nv).eq.2) then  !  Cloud base: conditional mean
                do i =1,iens
                 if(apoint(i).le.0.) apoint(i)=20000.0
                end do
                call get_cond_mean (apoint,iens,20000.0,amean,aspread,
     +                               miss,weight)

           else if(k4(nv).eq.3.and.k5(nv).eq.5.and.k6(nv).eq.3) then  !  Cloud top: conditional mean
                do i =1,iens
                 if(apoint(i).le.0.) apoint(i)=20000.0
                end do
                call get_cond_mean (apoint,iens,20000.0,amean,aspread,
     +                               miss,weight)

           else
             
             call getmean(apoint,iens,amean,aspread,miss,weight)

           end if

            vrbl_mn(igrid,lv)=amean
            vrbl_sp(igrid,lv)=aspread

            end do 
  
         !write(*,*)'Direct MEAN',(vrbl_mn(i,lv),i=20001,20010)
         !write(*,'(a4,10f9.2)')'SPRD',(vrbl_sp(i,lv),i=10001,10010)

           end do  !end if Mlvl

        END IF

        IF (trim(Psignal(nv)).eq.'P') THEN

         do lv = 1, Plvl(nv)
          do lt = 1, Tlvl(nv)

           do igrid=1,jf

            apoint = rawdata_pr(igrid,:,lv)
            miss = missing(nv,:)

          if ((k4(nv).eq.2.and.k5(nv).eq.222).or.
     +        (k4(nv).eq.2.and.k5(nv).eq.223).or.
     +        (k4(nv).eq.16.and.k5(nv).eq.198).or.
     +        (k4(nv).eq.16.and.k5(nv).eq.195).or.
     +        (k4(nv).eq.16.and.k5(nv).eq.196).or.
     +        (k4(nv).eq.7.and.k5(nv).eq.199).or.
     +        (k4(nv).eq.2.and.k5(nv).eq.220).or.
     +        (k4(nv).eq.2.and.k5(nv).eq.221) ) then
             do ie=1,iens
               if(mbrname(ie).eq.'hrrrgsd') then
                miss(ie)=1
               end if
             end do
          end if


            if(k4(nv).eq.19.and.k5(nv).eq.0) then  ! Visibility: conditional mean 
             do ie =1,iens
               if(apoint(ie).lt.1.0) miss(ie)=1    !exclude spc runs's error vis 
             end do
            end if

            if(k4(nv).eq.3.and.k5(nv).eq.5.and.k6(nv).eq.2) then  !  Cloud base: conditional mean
                do i =1,iens
                 if(apoint(i).le.0.) apoint(i)=20000.0
                end do
            end if
            if(k4(nv).eq.3.and.k5(nv).eq.5.and.k6(nv).eq.3) then  !  Cloud top: conditional mean
                do i =1,iens
                 if(apoint(i).le.0.) apoint(i)=20000.0
                end do
            end if


            if(trim(op(nv)).ne.'-') then

             thr1 = Thrs(nv,lt)
             thr2 = 0.
             call getprob(apoint,iens,thr1,thr2,op(nv),aprob,
     +                         miss,weight)
             vrbl_pr(igrid,lv,lt)=aprob

            else

             if(lt.lt.Tlvl(nv)) then
               thr1 = Thrs(nv,lt)
               thr2 = Thrs(nv,lt+1)
               call getprob(apoint,iens,thr1,thr2,op(nv),aprob,
     +                      miss,weight)
               vrbl_pr(igrid,lv,lt)=aprob
             end if

            end if

           end do

          if(vname(nv).eq.'AP3h') then    !This var will be used in cptp nightning/thunder computation
            if(thr1.eq.0.25) then
             p03mp01(:)=vrbl_pr(:,lv,lt)
            end if
          end if

c        write(*,'(a4,10f9.2)')'PROB',(vrbl_pr(i,lv,lt),i=10001,10010)

          end do !end of Tlvl
         end do !end of Plvl


        END IF


        write(*,*) 'Ensemble computation done for direct var ', nv
c Loop 1-3:  Packing  mean/spread/prob for this direct variable

       !ceiling/cloud base/freezing level can not use NAM's gfld, only can use RAP's !gfld
       !Other only can use RAP's gfld

       !On Dell VIS only can use its own gfld 

       if(k4(nv).eq.3.and.k5(nv).eq.5) then
          gfld_temp=gfld_ceil
       else if(k4(nv).eq.19.and.k5(nv).eq.0.and.k6(nv).eq.1) then
          gfld_temp=gfld_vis
       else
          gfld_temp=gfld_nam
       end if

      call packGB2_mean(imean,isprd,vrbl_mn,vrbl_sp,   !jpd12 is determined inside 
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid,gfld_temp)           !gfld is used to send in other info 

         write(*,*) 'packing mean direct var for', nv

      call packGB2_prob(iprob,vrbl_pr,             !jpd12 is determined inside
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lth,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid,gfld_temp)

         write(*,*) 'packing prob direct var for', nv

 222    CONTINUE 

c Loop 1-4: Deallocation

        if (allocated(rawdata_mn)) deallocate (rawdata_mn)
        if (allocated(rawdata_pr)) deallocate (rawdata_pr)
        if (allocated(vrbl_mn)) deallocate (vrbl_mn)
        if (allocated(vrbl_sp)) deallocate (vrbl_sp)
        if (allocated(vrbl_pr)) deallocate (vrbl_pr)

        write(*,*) 'Variable # ', nv, ' done ....'
      
2001     continue  !end of direct variables loop



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Loop 2: for derived variables
c    2-0: allocate necessary arrays
c    2-1: compute all derived variable by calling their correspondent subroutines 
c            Note: all input data will be read inside each subroutine 
c    2-2: pack derived varaibles
c    2-3: deallocation

       write(*,*) 'Now begin to generate derived products ......'

       DO 2002 nv = 1, nderiv

c   Loop 2-0: allocation
         gfld%discipline=0   !reset discipline# in case it was changed (i.e. fire weather)

         jpd1=dk4(nv)
         jpd2=dk5(nv)
         jpd10=dk6(nv)
         jpd27=-9999

          Lm=max(1,dMlvl(nv))
          Lp=max(1,dPlvl(nv))
          Lth=max(1,dTlvl(nv))
c          write(*,*) 'dTlvl(nv)=',dTlvl(nv),Lth

          if (.NOT.allocated(derv_mn)) then
            allocate (derv_mn(jf,Lm))         
          end if
          if (.NOT.allocated(derv_sp)) then
            allocate (derv_sp(jf,Lm))
          end if
          if (.NOT.allocated(derv_pr)) then
            allocate (derv_pr(jf,Lp,Lth))
          end if

c Loop 2-1: Computation 

cc%%%%%%% 1. To  see if there is thickness computation, if yes, do it
c           if (dk5(nv).eq.7.and. dk6(nv).eq.101) then
c

cc%%%%%%% 2. To see if there is precipitation type computation, if yes, do it
          if (dk4(nv).eq.1.and.dk5(nv).eq.19) then

              if (.NOT.allocated(ptype_mn)) then
                allocate (ptype_mn(jf,4))
              end if
              if (.NOT.allocated(ptype_pr)) then
                allocate (ptype_pr(jf,4))
              end if
            
              call preciptype (nv,ifunit,jf,iens,
     +         ptype_mn,ptype_pr)

              do jp=1,4
               jpd1=1
               jpd2=jptyp2(jp)
               jpd10=1
               jpd12=0
               jpd27=-999
               derv_mn(:,1)=ptype_mn(:,jp)     !for precip type, Lm,lp,Lth all are 1
               derv_sp=0.0
               derv_pr(:,1,1)=ptype_pr(:,jp)

               gfld_temp=gfld_nam
               call packGB2_mean_derv(imean,isprd,derv_mn,
     +              derv_sp,nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +              iens,iyr,imon,idy,ihr,ifhr,gribid,gfld_temp)  !gfld uses what was got from previous direct variables 

               gfld_temp=gfld_nam                           !some of idrtmpl() fields have been changed after packGB2_prob,
               call packGB2_prob_derv(iprob,derv_pr,
     +              nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lth,
     +              iens,iyr,imon,idy,ihr,ifhr,gribid,gfld_temp)  !gfld uses what was got from previous direct variables

              end do

              if (allocated(ptype_mn)) deallocate (ptype_mn)
              if (allocated(ptype_pr)) deallocate (ptype_pr)

c              write(*,*) 'preciptype done!'
            end if   !end of preciptype computation and save

cc%%%%%%% 3. To see if there is wind speed computation, if yes, do it
          if (dk4(nv).eq.2.and.dk5(nv).eq.1.and.dk6(nv).ne.108) then
            call wind (nv,ifunit,jf,iens,Lm,Lp,Lth,
     +        derv_mn,derv_sp,derv_pr,weight,mbrname)
c            write(*,*) 'Wind done'
          end if


cc%%%%%%% 4. To see if there is icing computation, if yes, do it
          if (dk4(nv).eq.19.and.dk5(nv).eq.7) then
            call get_icing (nv,ifunit,jf,iens,Lp,Lth,
     +              derv_pr,weight)
c            write(*,*) 'Icing done'
          end if

cc%%%%%%% 5. To see if there is CAT computation, if yes, do it
          if (dk4(nv).eq.19.and.dk5(nv).eq.22) then
            call get_cat (nv,ifunit,jf,iens,Lp,Lth,
     +              derv_pr,im,jm,dx,dy,weight)
            dPlvl(nv)=dPlvl(nv)-1
c            write(*,*) 'CAT done' 
          end if

cc%%%%%%% 6. To see if there is flight restriction  computation, if yes, do it
          if (dk4(nv).eq.19.and.dk5(nv).eq.205) then
            call  flight_res (nv,ifunit,jf,iens,Lp,Lth,
     +              derv_pr,weight)
c             write(*,*) 'Flight restiction done'
          end if

cc%%%%%%% 7. To see if there is Hains index for fire weather computation, if yes, do it
          if (dk4(nv).eq.4.and.dk5(nv).eq.2) then
            dTlvl(nv)=dTlvl(nv)-1
            Lth=Lth-1                !since dop='-'
            call fire_weather (nv,ifunit,jf,iens,Lp,Lth,
     +              derv_pr,weight)
            
            gfld%discipline=2      !Fireweather discipline = 2, used for packing
c            write(*,*) 'Hains Index done'
          end if

cc%%%%%%% 8. To see if there is ceiling computation, if yes, do it

          if (dk4(nv).eq.3.and.dk5(nv).eq.5.and.dk6(nv).eq.215) then
            call  getceiling (nv,ifunit,jf,iens,Lm,Lp,Lth,
     +             derv_mn,derv_sp,derv_pr,weight)
           
c            write(*,*) 'getceiling done', derv_mn(10000,1)
          end if


cc%%%%%%% 9. To see if there is fog  computation, if yes, do it

          if(dk4(nv).eq.6.and.dk5(nv).eq.193.and.dk6(nv).eq.103.
     +                             and.itime.ge.2) then

           call new_fog(nv,ifunit,ipunit,jf,im,jm,dx,dy,interval,
     +             iens,Lm,Lp,Lth,derv_mn,derv_sp,derv_pr,weight,
     +             ierr_open_punit)
c           write(*,*) 'new_fog done'

         end if

cc%%%%%%% 10. To see if there is thickness computation, if yes, do it
          if(dk4(nv).eq.3.and.dk5(nv).eq.5.and.
     +                             dk6(nv).eq.101) then

            call thickness (nv,ifunit,jf,iens,Lm,Lp,Lth,
     +             derv_mn,derv_sp,derv_pr,weight)

c              write(*,*) 'Thickness done'
          
          end if


cc%%%%%%% 11. To see if there is LLWS computation, if yes, do it
        if(dk4(nv).eq.2.and.dk5(nv).eq.192.and.
     +                             dk6(nv).eq.1) then


         !write(*,*) 'call llws ', eps       !??? Can not write to std !output   

         call llws (nv,ifunit,jf,iens,Lm,Lp,Lth,eps,
     +          derv_mn,derv_sp,derv_pr,weight)

c         write(*,*) 'LLWS done'

        end if

cc%%%%%%% 11. To see if there is CONVECTION  computation, if yes, do it
        if(dk4(nv).eq.1.and.dk5(nv).eq.196.and.
     +                             dk6(nv).eq.200) then
         call getconv(nv,ipunit,jf,im,jm,est,iens,Lm,Lp,Lth,
     +      gribid,derv_pr,weight,ierr_open_punit)

c           write(*,'(a9,10f9.2)')'CONV PROB',
c     +           (derv_pr(i,1,1),i=10001,10010)
c            write(*,*) 'getconv done' 

        end if


cc%%%%%%% 12. To see if there is Lighning, if yes, do it
c         missing() array is not considered, so current NSSE won't do this
         if(dk4(nv).eq.17.and.dk5(nv).eq.192.and.
     +                             dk6(nv).eq.1) then
         call get_cptp_severe(nv,cycle(ihr+1),cfhr,ifunit,p03mp01,
     +       im,jm,km,jf,iens,Lp,Lth,derv_pr,weight,eps)
c           write(*,*) 'CPTP done '
         end if


cc%%%%%%% 13. To see if there is 850-300mb mean wind, if yes, do it
         if(dk4(nv).eq.2.and.dk5(nv).eq.1.and.
     +                             dk6(nv).eq.108) then

           call meanwind(nv,ifunit,jf,iens,Lm,Lp,Lth,eps,
     +          derv_pr,weight,mbrname)

c            write(*,*) '850-300mb mean-wind done'
         end if




C  Loop 2-2:  Pack Non-precip type dereved mean/spread/prob products
C
          if(dk4(nv).eq.1.and.dk5(nv).eq.19 ) then  !Precip type prob is packed already
             goto 20021
          else                                       !Non-precip type
               jpd1=dk4(nv)
               jpd2=dk5(nv)
               jpd10=dk6(nv)
               jpd27=-9999
               !write(*,*) 'pack _mean_derv for ', nv

               !Ceiling use RAP's gfld, other use RAP's gfld
               if(dk4(nv).eq.3.and.dk5(nv).eq.5)then
                gfld_temp=gfld_ceil
               else
                gfld_temp=gfld_nam
               end if

               write(*,*) 'pack _derv for ', nv
     +          ,jpd1,jpd2,jpd10,jpd27,iprob,jf,Lp,Lth,
     +          iens,iyr,imon,idy,ihr,ifhr,gribid

               call packGB2_prob_derv(iprob,derv_pr,
     +              nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lth,
     +              iens,iyr,imon,idy,ihr,ifhr,gribid,gfld_nam)  !borrow gfld from what was got from previous direct variables

               call packGB2_mean_derv(imean,isprd,derv_mn,
     +              derv_sp,nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +              iens,iyr,imon,idy,ihr,ifhr,gribid,gfld_temp) !borrow gfld from what was got from previous direct variables

               write(*,*) 'pack _derv done for', nv


          end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Loop 2-3: deallocation

20021    continue

       if (allocated(derv_mn)) deallocate (derv_mn)
       if (allocated(derv_sp)) deallocate (derv_sp)
       if (allocated(derv_pr)) deallocate (derv_pr)




2002   CONTINUE  ! end of derived variables loop


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c Loop 3: Process MXP variables (min,max, 10%,20%, ...90%) 
c 3-0: Allocate array mxp8 
c 3-1: Compute each mxp variable 
c 3-3: pack all mxp variables
c 3-4: deallocate
     
c       write(*,*) 'Now begin to generate MXP products ......'
       
       if(nmxp.gt.0) then
         pmin=trim(eps)//'.pmin.t'//cycle(ihr+1)//'z'//'.f'//cfhr
         pmax=trim(eps)//'.pmax.t'//cycle(ihr+1)//'z'//'.f'//cfhr
         pmod=trim(eps)//'.pmod.t'//cycle(ihr+1)//'z'//'.f'//cfhr
         pp10=trim(eps)//'.pp10.t'//cycle(ihr+1)//'z'//'.f'//cfhr
         pp25=trim(eps)//'.pp25.t'//cycle(ihr+1)//'z'//'.f'//cfhr
         pp50=trim(eps)//'.pp50.t'//cycle(ihr+1)//'z'//'.f'//cfhr
         pp75=trim(eps)//'.pp75.t'//cycle(ihr+1)//'z'//'.f'//cfhr
         pp90=trim(eps)//'.pp90.t'//cycle(ihr+1)//'z'//'.f'//cfhr

         call baopen(301,pmin,ierr)
         call baopen(302,pmax,ierr)
         call baopen(303,pmod,ierr)
         call baopen(304,pp10,ierr)
         call baopen(305,pp25,ierr)
         call baopen(306,pp50,ierr)
         call baopen(307,pp75,ierr)
         call baopen(308,pp90,ierr)
       end if

       DO 2003 nv = 1, nmxp

c  Loop 3-0: allocation 

         Lq=max(1,qMlvl(nv))
         if (.NOT.allocated(mxp8)) then
          allocate (mxp8(jf,Lq,8))
         end if
c Loop 3-1: Computation
         call get_mxp(nv,ifunit,jf,iens,Lq,mxp8,weight) 
         write(*,*) 'call get_mxp done for ', nv

c Loop 3-2: Packing 

         gfld_temp=gfld_nam

         jpd1=qk4(nv)
         jpd2=qk5(nv)
         jpd10=qk6(nv)
         jpd27=-9999


         do kq=1,8
          call packGB2_mxp(iqout(kq),kq,mxp8,    
     +      nv,jpd1,jpd2,jpd10,jpd27,jf,Lq,
     +      iens,iyr,imon,idy,ihr,ifhr,gribid,gfld_temp)           !gfld_temp  is used to send in other info 
         end do


        if (allocated(mxp8)) deallocate(mxp8)
2003  CONTINUE
       

       do irun=1,iens
        call baclose(ifunit(irun),ierr)
        call baclose(ipunit(irun),ierr)
       end do

       call baclose(imean,ierr)
       call baclose(isprd,ierr)
       call baclose(iprob,ierr)

       if (nmxp.gt.0) then
         call baclose(301,ierr)
         call baclose(302,ierr)
         call baclose(303,ierr)
         call baclose(304,ierr)
         call baclose(305,ierr)
         call baclose(306,ierr)
         call baclose(307,ierr)
         call baclose(308,ierr)
       end if
 
3000    continue   ! end of itime loop
      

      write(*,*) 'Entire program completed!'

      stop
      end



