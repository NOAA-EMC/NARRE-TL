

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine getceiling: compute ceiling height
c     
c     Author: Binbin Zhou, Jun 1, 2006
c     Modification: 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine getceiling (nv,ifunit,jf,iens,Lm,Lp,Lt,  
     +             derv_mn,derv_sp,derv_pr,wgt)
        
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
        REAL,dimension(jf,Lm),intent(INOUT) :: derv_mn
        REAL,dimension(jf,Lm),intent(INOUT) :: derv_sp
        REAL,dimension(jf,Lp,Lt),intent(INOUT) :: derv_pr

        REAL, dimension(jf,iens) :: tcld,cldb,hsfc

        real CEILapoint(iens),TCLDapoint(iens),CLDBapoint(iens),
     +       HSFCapoint(iens)
        real wgt(30), amean,aspread,x(iens)

        integer miss(iens)

        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld


        jpdtn=0
        jp27=-9999

        write(*,*) 'In getceiling .....'
        write(*,*) 'nv,ifunit,jf,iens,Lm,Lp,Lt',
     +              nv,ifunit,jf,iens,Lm,Lp,Lt

        miss=0
        do 400 irun=1,iens

           call readGB2(ifunit(irun),jpdtn,6,1,200,0,jp27,gfld,iret) !total cloud
            if(iret.eq.0) then
             tcld(:,irun)=gfld%fld
            else
             write(*,*) 'Try jpd10=10 for RAP' 
             call readGB2(ifunit(irun),jpdtn,6,1,10,0,jp27,gfld,iret) 
             if(iret.eq.0) then
              write(*,*) 'Read RAP total ok'
              tcld(:,irun)=gfld%fld
             end if
            end if

           call readGB2(ifunit(irun),jpdtn,3,5,2,0,jp27,gfld,iret)   !cloud base
            if(iret.eq.0) then            
             cldb(:,irun)=gfld%fld
            else
             miss(irun)=1
             goto 400
            end if

           call readGB2(ifunit(irun),jpdtn,3,5,1,0,jp27,gfld,iret)  !sfc height 
            if(iret.eq.0) then
             hsfc(:,irun)=gfld%fld
            else
             miss(irun)=1
             goto 400
            end if

 400    continue

        do 600 igrid = 1,jf

          TCLDapoint=tcld(igrid,:)
          CLDBapoint=cldb(igrid,:)
          HSFCapoint=hsfc(igrid,:)

          do i = 1, iens
            if(miss(i).eq.0) then
              if(CLDBapoint(i).lt.0.0) CEILapoint(i)=20000.0    !Dec. 30, 'le'->'lt'
              if(TCLDapoint(i).ge.50.0 .and.
     +          CLDBapoint(i).ge.0.0  ) then
                CEILapoint(i) = CLDBapoint(i) - HSFCapoint(i)
                if(CEILapoint(i).lt.0.0) CEILapoint(i)=0.0
              else
                CEILapoint(i)=20000.0
              end if
            end if
          end do               

          call get_cond_mean (CEILapoint, iens, 
     +      20000.0, amean, aspread,miss,wgt)
           derv_mn(igrid,1)=amean
           derv_sp(igrid,1)=aspread

            do 30 lv=1,dPlvl(nv)
              do lt = 1, dTlvl(nv)

                if(trim(dop(nv)).ne.'-') then
                 thr1 = dThrs(nv,lt)
                 thr2 = 0.
                 call getprob(CEILapoint,iens,
     +                thr1,thr2,dop(nv),aprob,miss,wgt)
                 derv_pr(igrid,lv,lt)=aprob
                else
                  if(lt.lt.dTlvl(nv)) then
                    thr1 = dThrs(nv,lt)
                    thr2 = dThrs(nv,lt+1)
                    call getprob(CEILapoint,iens,
     +                   thr1,thr2,dop(nv),aprob,miss,wgt)
                    derv_pr(igrid,lv,lt)=aprob
                   end if
                end if
             end do

30         continue  

600      continue

        return
        end

       subroutine get_cond_mean_ceil (x,n,alarge,mean,spread,miss,wgt)
         real x(n),wgt(30), mean, spread, alarge
         integer n,miss(n)

         mean = 0.
         count = 0.

         wsum=0.0
         count = 0.
         do 102 i=1,n
           if(miss(i).eq.0) then
             if(x(i).ge.alarge) goto 102
             wsum=wsum+wgt(i)
             count = count + 1.0
           end if
102      continue

         do 100 i=1,n
          if (miss(i).eq.0) then
           if(x(i).ge.alarge) goto 100
            mean = mean + x(i)*(wgt(i)/wsum)
          end if
100      continue
         if(count.gt.0.0) then
          mean=mean/count
         else
          mean=alarge
         end if

         spread = 0.
         do 200 i = 1, n
          if(x(i).ge.alarge.or.miss(i).eq.1) goto 200                       !!!  modified in Apr. 9 2008 
           spread = spread + (wgt(i)/wsum)*(x(i)-mean)**2
200      continue

         if (mean.eq.alarge) then
             spread=0.0
         else
             spread = sqrt (spread )
         end if

         return
         end

