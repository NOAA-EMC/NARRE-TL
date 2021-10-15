cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine getfog: compute fog occurrence
c     
c     Author: Binbin Zhou, Aug 1, 2009
c     Modification: adapted from Beijing Olympic Project's product generator
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine getfog (nv,rawdata,jf,iens,
     +                     derv_pr,missing,wgt)

            include 'parm.inc'

C for variable table:
        Integer numvar
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
        Integer dk5(maxvar), dk6(maxvar), dk4(maxvar)
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
                                                                                                                                                                

        INTEGER, intent(IN) :: nv, jf, iens
        REAL,dimension(jf,iens,maxmlvl),intent(IN) :: rawdata
        REAL,dimension(jf,maxplvl,maxtlvl),intent(INOUT) :: 
     +               derv_pr

        real FOGapnt(iens),CLDBapnt(iens),CLDTapnt(iens),
     +       U10apnt(iens),V10apnt(iens),RH2apnt(iens),W10apnt(iens),
     +       VISapnt(iens),SFCapnt(iens)

        real fog(jf,iens),aprob, lwc(jf,iens), amean
        real wgt(30)

        integer,dimension(jf,iens),intent(IN) :: missing 
        integer miss(iens)      


        ID_CLDB = index_table(k5,k6,7,2,maxvar)      !search index of cloud base heightin the table
        ID_CLDT = index_table(k5,k6,7,3,maxvar)      !search index of cloud cloud heightin the table
        ID_U10   = index_table(k5,k6,33,105,maxvar)  !search index of 10m U
        ID_V10   = index_table(k5,k6,34,105,maxvar)  !search index of 10m V
        ID_RH2   = index_table(k5,k6,52,105,maxvar)  !search index of 2m RH
        ID_VIS   = index_table(k5,k6,20,1,maxvar)    !search index of sfc visibility       
        ID_SFC   = index_table(k5,k6,7,1,maxvar)     !search index of sfc height in the table

        write(*,*) 'In getfog -->'

        write(*,*)'ID_U10, ID_V10, ID_RH2, ID_VIS=',
     +             ID_U10, ID_V10, ID_RH2,ID_VIS
        write(*,*)'ID_CLDB, ID_CLDT, ID_SFC',
     +    ID_CLDB, ID_CLDT, ID_SFC

      if (ID_CLDB .gt.0 .and. ID_CLDT .gt.0 .and. 
     +    ID_U10  .gt.0 .and. ID_V10  .gt.0 .and. 
     +    ID_RH2  .gt.0 .and. ID_VIS  .gt.0 .and. 
     +    ID_SFC  .gt.0 ) then

        do igrid = 1,jf

          CLDBapnt=rawdata(igrid,:,ID_CLDB,1)
          CLDTapnt=rawdata(igrid,:,ID_CLDT,1)
          RH2apnt =rawdata(igrid,:,ID_RH2,1)
          W10apnt=sqrt(rawdata(igrid,:,ID_U10,1)**2 +
     +                 rawdata(igrid,:,ID_V10,1)**2 )         
          VISapnt =rawdata(igrid,:,ID_VIS,1)
          SFCapnt =rawdata(igrid,:,ID_SFC,1)

          if(igrid.eq.739400) then
            write(*,*) CLDBapnt
            write(*,*) CLDTapnt
            write(*,*) W10apnt
            write(*,*) RH2apnt 
            write(*,*) VISapnt 
            write(*,*) SFCapnt 
          end if

          FOGapnt = 0.
          hasfog= 0.

          do i = 1, iens

            if(CLDBapnt(i).lt.0.0) CLDBapnt(i)=20000.0
            if(CLDTapnt(i).lt.0.0) CLDTapnt(i)=20000.0

            if(CLDBapnt(i).ge.0.0) Then
              CLDBapnt(i) = CLDBapnt(i) - SFCapnt(i)
              if(CLDBapnt(i).lt.0.0)  CLDBapnt(i) = 0.0
            end if

            if(CLDTapnt(i).ge.0.0) Then
             CLDTapnt(i) = CLDTapnt(i) - SFCapnt(i)
             if(CLDTapnt(i).le.0.0) CLDTapnt(i) = 0.0
            end if

           if(CLDTapnt(i).LE.CLDBapnt(i)) CLDTapnt(i)=CLDBapnt(i)

             if(((CLDBapnt(i).ge.0.0 .and. CLDBapnt(i).le.50.0).and.  !50:  vertical resolution ~ 30m and detect "foot fog"
     +         (CLDTapnt(i).ge.0.0 .and. CLDTapnt(i).le.400.0))       !400: detect sea fog, advection fog and deep fogs
     +       .or. (W10apnt(i).le.0.05 .and. RH2apnt(i).ge.99.0)        !RH-wind to detect ground fog
     +       .or. VISapnt(i).le.1000.  ) then                         !equavalent to LWC = 0.015 g/kg

                FOGapnt(i) = 1.
                hasfog=1.

                if(RH2apnt(i).lt.95.0) FOGapnt(i)=0.                     !if RH<95%, no fog anyway
             else
                FOGapnt(i) = 0.
             endif

          end do

            fog(igrid,:)=FOGapnt(:)

        end do

           do 30 lv=1,dPlvl(nv)

             derv_pr(:,nv,lv,:)=0.0

             do lt = 1, dTlvl(nv)
 
              thr1 = dThrs(nv,lt)
                             
              do igrid = 1,jf

               miss=missing(igrid,:)

               FOGapnt=fog(igrid,:)

               call getprob(FOGapnt,iens,
     +                thr1,thr2,dop(nv),aprob,miss,wgt)
          
               derv_pr(igrid,nv,lv,lt)=aprob 
        
              end do
             end do

30          continue  

        end if

        return
        end
