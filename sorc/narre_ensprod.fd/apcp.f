cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Subroutine apcp: compute 3/6/12/24 hour accumulated precipitation
c  mean/sprea/probability based on Jun Du's old version
c  Note: The mean is "unditional mean", ie. 0 precipitation is also
c        taken into account in mean
c
c  Author: Binbin Zhou, Agu. 8, 2005
c
c  Input: nv, itime, i00, precip, jf, iens, interval
c  output: derv_mn,  derv_pr
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   	subroutine apcp(nv,itime,i00,precip,jf,iens,interval,
     +   loutput, derv_mn,derv_sp,derv_pr,prcpmax,prcpmin,
     +   missing,wgt) 

         include 'parm.inc'

C for variable table:
        Integer numvar, nderiv
        Character*4 vname(maxvar)
        Integer k5(maxvar), k6(maxvar),k4(maxvar)
        Character*1 Msignal(maxvar), Psignal(maxvar)
        Integer Mlvl(maxvar), MeanLevel(maxvar,maxmlvl)
        Integer Plvl(maxvar), ProbLevel(maxvar,maxplvl)
        Integer Tlvl(maxvar)
        Character*1 op(maxvar)
        Real    Thrs(maxvar,maxtlvl)
                                                                                                                                                                
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
                                                                                                                                                                
        common /tbl/numvar,
     +              vname,k4,k5,k6,Mlvl,Plvl,Tlvl,
     +              MeanLevel,ProbLevel,Thrs,
     +              Msignal,Psignal,op
                                                                                                                                                                
        common /dtbl/nderiv,
     +              dvname,dk4,dk5,dk6,dMlvl,dPlvl,dTlvl,
     +              dMeanLevel,dProbLevel,dThrs,
     +              dMsignal,dPsignal,MPairLevel,PPairLevel,dop


        INTEGER, intent(IN) :: nv,itime,i00,jf,iens,interval,loutput
        REAL,dimension(jf,iens,loutput),intent(IN) :: precip
        REAL,dimension(jf,maxmlvl),intent(INOUT) :: derv_mn
        REAL,dimension(jf,maxmlvl),intent(INOUT) :: derv_sp
        REAL,dimension(jf,maxplvl,maxtlvl),intent(INOUT) :: 
     +               derv_pr
        REAL,dimension(jf,maxmlvl),intent(INOUT)        :: prcpmax
        REAL,dimension(jf,maxmlvl),intent(INOUT)        :: prcpmin
       

        INTEGER               :: tick ! counter for counting previous times           
        REAL,dimension(iens)  :: apoint
        REAL                  :: amean,aspread,aprob
 
        REAL wgt(30),miss(iens)
        INTEGER,dimension(jf,iens),intent(IN) :: missing
         
          ID_PRCP = index_table(k5,k6,dk5(nv),dk6(nv),maxvar)       !search index of direct variable in the table
 
 
          if (ID_PRCP .gt.0 ) then 

             write(*,*) 'In APCP-->',dMlvl(nv),dPlvl(nv)

             do lv=1,dMlvl(nv)                                  !for accumulated hours
               write(*,*) itime, dMeanLevel(nv,lv)      
               if(itime.ge.dMeanLevel(nv,lv)) then

                 do igrid = 1,jf
                   apoint = 0.
                   miss=missing(igrid,:)

                   tick = dMeanLevel(nv,lv)/interval - 1
                   do while ( tick .ge. 0 )
                     apoint= apoint + precip(igrid,:,i00-tick)
                     tick = tick - 1
                   end do
                   call getmean(apoint,iens,amean,aspread,miss,wgt)
                   derv_mn(igrid,nv,lv) = amean
                   derv_sp(igrid,nv,lv) = aspread
                   prcpmax(igrid,lv) = maxval(apoint)
                   prcpmin(igrid,lv) = minval(apoint)

                 end do
               end if
             end do

             do lv=1,dPlvl(nv)

              write(*,*) itime, dProbLevel(nv,lv)
              if(itime.ge.dProbLevel(nv,lv)) then

                do lt =1,dTlvl(nv)
                  do igrid = 1, jf
                    apoint = 0.
                    miss=missing(igrid,:)

                    tick = dProbLevel(nv,lv)/interval - 1
                    do while ( tick .ge. 0 )
                     apoint= apoint + precip(igrid,:,i00-tick)
                     tick = tick - 1
                    end do

                    if(trim(dop(nv)).ne.'-') then
                     thr1 = dThrs(nv,lt)
                     thr2 = 0.
               call getprob(apoint,iens,thr1,thr2,dop(nv),aprob,
     +               miss,wgt)
                     derv_pr(igrid,nv,lv,lt)=aprob
                    else
                     if(lt.lt.dTlvl(nv)) then
                       thr1 = dThrs(nv,lt)
                       thr2 = dThrs(nv,lt+1)
               call getprob(apoint,iens,thr1,thr2,dop(nv),aprob,
     +               miss,wgt)
                       derv_pr(igrid,nv,lv,lt)=aprob
                     end if
                    end if

                  end do
                end do

              end if
             end do
           
           end if

           return
           end
