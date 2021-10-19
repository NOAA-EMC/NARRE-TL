cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine thickness: compute distance between 2 pressure levels
c     first search for index of HGT from the table -> IDV, if found:
c     then search level pair (MPairlevel) index from MeanLevel-> L1, L2, if found:
c     then compute the difference between the two levels
c     
c     Author: Binbin Zhou, Aug, 5, 2005
c     05/15/2013: Updated to grib2 I/O, B. Zhou
c
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine thickness (nv,ifunit,jf,iens,Lm,Lp,Lt,  
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
        REAL,dimension(jf,Lm,Lt),intent(INOUT) :: derv_pr

        INTEGER miss(iens)

        real apoint(iens),wgt(30)
        REAL, dimension(jf,iens) :: h1,h2
        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld


        jpdtn=0
        jdp1=dk4(nv)
        jdp2=dk5(nv)
        jpd10=100       
        jp27=-9999


c        write(*,*) 'In thickness .....'
c        write(*,*) 'nv,ifunit,jf,iens,Lm,Lp,Lt,jpd10',
c     +              nv,ifunit,jf,iens,Lm,Lp,Lt,jpd10

         miss=0
          do 400 irun=1,iens
           jpd12=MPairLevel(nv,1,1)
           call readGB2(ifunit(irun),jpdtn,jdp1,jdp2,jpd10,jpd12,
     +                                          jp27,gfld,ie) !higher level  
           if(ie.eq.0) then
             h1(:,irun)=gfld%fld
           else
             miss=1
             goto 400
           end if

           jpd12=MPairLevel(nv,1,2)
           call readGB2(ifunit(irun),jpdtn,jdp1,jdp2,jpd10,jpd12,
     +                                          jp27,gfld,ie) !lower level  
           if(ie.eq.0) then
            h2(:,irun)=gfld%fld
           else
             miss=1
             goto 400
           end if

 400      continue



             do lv=1,dMlvl(nv)                       !for all pairs of levels

                do igrid = 1,jf

                  apoint = abs( h1(igrid,:) - h2(igrid,:))
                  call getmean(apoint,iens,amean,aspread,
     +                 miss,weight)
                  derv_mn(igrid,lv)=amean
                  derv_sp(igrid,lv)=aspread
                end do


              end do

              do lv=1,dPlvl(nv)
                do lt = 1, dTlvl(nv)

                 do igrid = 1,jf
                   apoint = abs( h1(igrid,:) - h2(igrid,:))

                    if(trim(dop(nv)).ne.'-') then
                     thr1 = dThrs(nv,lt)
                     thr2 = 0.
               call getprob(apoint,iens,thr1,thr2,dop(nv),aprob,
     +                       miss,weight)
                     derv_pr(igrid,lv,lt)=aprob
                    else
                     if(lt.lt.dTlvl(nv)) then
                       thr1 = dThrs(nv,lt)
                       thr2 = dThrs(nv,lt+1)
               call getprob(apoint,iens,thr1,thr2,dop(nv),aprob)
                       derv_pr(igrid,lv,lt)=aprob
                     end if
                    end if

                 end do

                end do
              end do


           return
           end
