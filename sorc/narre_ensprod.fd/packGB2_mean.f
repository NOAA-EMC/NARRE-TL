      subroutine packGB2_mean(imean,isprd,vrbl_mn,
     +     vrbl_sp,nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid,gfld)

        use grib_mod
        include 'parm.inc'

        type(gribfield) :: gfld

C for variable table:
        Integer numvar
        Character*4 vname(maxvar)
        Integer k5(maxvar), k6(maxvar), k4(maxvar)
        Character*1 Msignal(maxvar), Psignal(maxvar)
        Integer Mlvl(maxvar), MeanLevel(maxvar,maxmlvl)
        Integer Plvl(maxvar), ProbLevel(maxvar,maxplvl)
        Integer Tlvl(maxvar)
        Character*1 op(maxvar)
        Real    Thrs(maxvar,maxtlvl)


        common /tbl/numvar,
     +              vname,k4,k5,k6,Mlvl,Plvl,Tlvl,
     +              MeanLevel,ProbLevel,Thrs,
     +              Msignal,Psignal,op


         INTEGER,intent(IN) :: imean,isprd,
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid

        REAL,dimension(jf,Lm),intent(IN) :: vrbl_mn
        REAL,dimension(jf,Lm),intent(IN) :: vrbl_sp

        INTEGER,allocatable,dimension(:) ::   ipdtmpl

        integer ml

        !write(*,*) 'packing direct mean/spred for nv ', nv

c        write(*,*) nv, imean,isprd,
c     +     jpd1,jpd2,jpd10,jpd27,jf,Lm,
c     +     iens,iyr,imon,idy,ihr,ifhr,gribid



        DO 1000 ml=1,Mlvl(nv)

            !redefine some of gfld%idsect() array elements 
            gfld%idsect(1)=7
            gfld%idsect(2)=2
            gfld%idsect(3)=0   !experimental, see Table 1.0
            gfld%idsect(4)=1   !experimental, see Table 1.1
            gfld%idsect(5)=1   !experimental, see Table 1.2
            gfld%idsect(6)=iyr  !year
            gfld%idsect(7)=imon !mon
            gfld%idsect(8)=idy  !day
            gfld%idsect(9)=ihr  !cycle time
            gfld%idsect(10)=0
            gfld%idsect(11)=0
            gfld%idsect(12)=0
            gfld%idsect(13)=1


        !Pack direct variable mean/spread

         if (ml.eq.1) then              !this is just do once

            if (jpd1.eq.1.and.(jpd2.eq.8.or.jpd2.eq.11) ) then  
             ipdtnum=12                  !ensemble APCP mean use Template 4.12
             ipdtlen=31
            else 
             ipdtnum=2                   !ensemble NON-accum mean use Template 4.2
             ipdtlen=17
            end if

            allocate (ipdtmpl(ipdtlen))

             ipdtmpl(1)=jpd1
             ipdtmpl(2)=jpd2
             ipdtmpl(3)=4
             ipdtmpl(4)=0
             ipdtmpl(5)=117              !shared with NARRE-TL suggested by Tom Hultquist
             ipdtmpl(6)=0
             ipdtmpl(7)=0
             ipdtmpl(8)=1
             ipdtmpl(9)=ifhr             !fcst time    
             ipdtmpl(10)=jpd10
             ipdtmpl(11)=gfld%ipdtmpl(11)
             !ipdtmpl(12)= see below 
             ipdtmpl(13)=gfld%ipdtmpl(13)
             ipdtmpl(14)=gfld%ipdtmpl(14)
             ipdtmpl(15)=gfld%ipdtmpl(15)
             ipdtmpl(16)=1               !weighted mean of all members
             ipdtmpl(17)=iens            !number of members
         
            if (jpd1.eq.1.and.(jpd2.eq.8.or.jpd2.eq.11) ) then  !Template 4.12 has extra elements than Template 4.2
              !2015121205 correction: B. zhou ...
              !ihr_ifhr=ihr+ifhr
              !call get_ymd(iyr,imon,idy,ihr_ifhr,kyr,kmon,kdy,khr)
              !ipdtmpl(18)=iyr   !year
              !ipdtmpl(19)=imon  !mon
              !ipdtmpl(20)=idy   !day
              !ipdtmpl(9) =ihr+ifhr-jpd27    !overwrite for APCP:
              !Beginning time of accumulation
              !ipdtmpl(21)=ihr+ifhr         !end time of accumulation
              !ipdtmpl(21)=khr         !end time of accumulation

              call get_time_GB2(iyr,imon,idy,ihr,ifhr,
     +            iyr1,imon1,idy1,ihr1)             
              ipdtmpl(18)=iyr1   !year
              ipdtmpl(19)=imon1  !mon
              ipdtmpl(20)=idy1   !day
              ipdtmpl(9)= ifhr-jpd27 
              ipdtmpl(21)=ihr1

              ipdtmpl(22)=0   
              ipdtmpl(23)=0  
              ipdtmpl(24)=1   
              ipdtmpl(25)=0   
              ipdtmpl(26)=1   
              ipdtmpl(27)=2                 !See Table 4.11, same start time (or same cycle time) 
              ipdtmpl(28)=1  
              ipdtmpl(29)=jpd27 
              ipdtmpl(30)=1
              ipdtmpl(31)=0   
             end if
 
          end if  

          if(jpd10.eq.100) then
             ipdtmpl(12)=MeanLevel(nv,ml)*100
          else
             ipdtmpl(12)=MeanLevel(nv,ml)
          end if

          gfld%fld=vrbl_mn(:,ml) 


          !!!IMPORTANT 
          !if(gfld%ibmap .ne. 255) then   ! exclude SREF mean 
          ! gfld%ibmap=0
          !end if
          if (gfld%idrtmpl(7).eq.255) gfld%idrtmpl(7)=0

          iret=0
          !write(*,*) 'ipdtnum, ipdtlen=',ipdtnum, ipdtlen
          !do k=1,ipdtlen
          !  write(*,*) k, ipdtmpl(k)
          !end do

          call Zputgb2(imean,gfld,ipdtmpl,ipdtnum,ipdtlen,iret)
          !if(iret.ne.0) then
          ! write(*,*) 'Zputgb2 mean error:',iret
          !end if

      
          gfld%fld=vrbl_sp(:,ml)
          ipdtmpl(16)=4                !Spread 

          iret=0
          call Zputgb2(isprd,gfld,ipdtmpl,ipdtnum,ipdtlen,iret)
          
          !if(iret.ne.0) then
          ! write(*,*) 'Zputgb2 spread error:',iret
          !end if
          
1000    CONTINUE   


        return
        end
