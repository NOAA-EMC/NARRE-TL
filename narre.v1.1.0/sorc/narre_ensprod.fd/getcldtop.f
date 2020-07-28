cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine getfog: compute mean cloud top height above the surface
c                        and probability
c   
c     Author: Binbin Zhou, Jun 1, 2006
c     Modification: 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine getcldtop (nv,rawdata,jf,iens,
     +      vrbl_mn, vrbl_sp, vrbl_pr,missing,wgt)

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
                                            
 
        common /tbl/numvar,
     +              vname,k4,k5,k6,Mlvl,Plvl,Tlvl,
     +              MeanLevel,ProbLevel,Thrs,
     +              Msignal,Psignal,op
                                                                                                                                                                                                                                                                        

        INTEGER, intent(IN) :: nv, jf, iens
        REAL,dimension(jf,iens,maxmlvl),intent(IN) :: rawdata
        REAL,dimension(jf,maxmlvl),intent(INOUT) :: vrbl_mn
        REAL,dimension(jf,maxmlvl),intent(INOUT) :: vrbl_sp
        REAL,dimension(jf,maxplvl,maxtlvl),intent(INOUT) ::
     +               vrbl_pr


        real CTOPapnt(iens),CLDTapnt(iens),SFCHapnt(iens)
        real ctop(jf,iens),aprob

        integer, dimension(jf,iens),intent(IN) :: missing
        integer miss(iens)


        ID_CLDT = index_table(k5,k6,7,3,maxvar)      !search index of cloud top height in the table
        ID_SFCH = index_table(k5,k6,7,1,maxvar)      !search index of surface height in the table

        write(*,*) 'In getcldtop -->'
      write(*,*)'ID_CLDT,ID_SFCH,Jf=',ID_CLDT,ID_SFCH,jf


      if (ID_CLDT .gt.0 .and. ID_SFCH .gt.0 ) then


        do igrid = 1,jf

           CLDTapnt=rawdata(igrid,:,ID_CLDT,1)
           SFCHapnt=rawdata(igrid,:,ID_SFCH,1)

           if(igrid.eq.2051) then
             write(*,*) 'CLDTapnt=',CLDTapnt
             write(*,*) 'SFCHapnt=',SFCHapnt
           end if


           do i = 1, iens
            if(CLDTapnt(i).ge.0.0  ) then                     !Dec. 30, 2008: 'gt'->'ge'
              CTOPapnt(i) = CLDTapnt(i) - SFCHapnt(i)
              if(CTOPapnt(i).lt.0.0) CTOPapnt(i)=0.0
            else
              CTOPapnt(i)=20000.0
            end if
           end do

          ctop(igrid,:)=CTOPapnt(:)

          miss=missing(igrid,:)

          call get_cond_mean (CTOPapnt,iens, 20000.0,
     +           amean, aspread,miss,wgt)

          vrbl_mn(igrid,nv,1)=amean
          vrbl_sp(igrid,nv,1)=aspread

        end do

        do 30 lv=1,Plvl(nv)

          do lt = 1, Tlvl(nv)
           do igrid = 1,jf

             CTOPapnt=ctop(igrid,:)

             miss=missing(igrid,:)

             if(trim(op(nv)).ne.'-') then
                 thr1 = Thrs(nv,lt)
                 thr2 = 0.
                 call getprob(CTOPapnt,iens,
     +         thr1,thr2,op(nv),aprob,miss,wgt)
                 vrbl_pr(igrid,nv,lv,lt)=aprob
              else
                if(lt.lt.Tlvl(nv)) then
                  thr1 = Thrs(nv,lt)
                  thr2 = Thrs(nv,lt+1)
                  call getprob(CTOPapnt,iens,
     +          thr1,thr2,op(nv),aprob,miss,wgt)
                  vrbl_pr(igrid,nv,lv,lt)=aprob
                end if
              end if

            end do
          end do

30          continue

      end if

      return
      end
