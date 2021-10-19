cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine dtra: compute variances from global ensemble system
c     first search for index of U, V from the table -> ID_U, ID_V,  
c         if found: then search level index from MeanLevel-> L, if found:
c     then compute variances as defined in the subroutine
c
c     Author: Binbin Zhou, Aug, 30, 2006
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine dtra (nv,rawdata, jf, iens,
     +             derv_dtra)
                                                                                                                                    
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
        Integer ID_U,ID_V,ID_W,ID_T,ID_Q
        Integer L1,L2,L3,L4,L5
 
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
        REAL,dimension(jf,maxmlvl,10),intent(INOUT) :: derv_dtra
 
        real apoint(iens),Uapoint(iens),Vapoint(iens),
     +    Wapoint(iens),Tapoint(iens),Qapoint(iens)

          if(dk6(nv).eq.105) then
           ID_U = index_table(k5,k6,33,105,maxvar)      !search index of direct variable U in the table
           ID_V = index_table(k5,k6,34,105,maxvar)      !search index of direct variable V in the table
          end if
 
          if(dk6(nv).eq.100) then
           ID_U = index_table(k5,k6,33,100,maxvar)      !search index of direct variable U in the table
           ID_V = index_table(k5,k6,34,100,maxvar)      !search index of direct variable V in the table
          end if
 
          if(dk6(nv).eq.109) then                       !for sigma level
           ID_U = index_table(k5,k6,33,109,maxvar)      !search index of direct variable U in the table
           ID_V = index_table(k5,k6,34,109,maxvar)      !search index of direct variable V in the table
          end if

           write(*,*) 'ID_U =',ID_U 
           write(*,*) 'ID_V =',ID_V 

          if (ID_U .gt.0 .and. ID_V .gt.0 ) then

             do lv=1,dMlvl(nv)                       !for all pairs of levels
 
                L1=index_int_array(MeanLevel(ID_U,:),           !search index of 1st level from
     +                    dMeanlevel(nv,lv),maxmlvl)            !   ID_U-th variable's MeanLevel
                if (L1.eq.0) then
                 write(*,*)dMeanlevel(nv,lv),' not found for U'
                 stop
                end if
 
                L2=index_int_array(MeanLevel(ID_V,:),           !search index of 1st level from
     +                    dMeanlevel(nv,lv),maxmlvl)            !   ID_U-th variable's MeanLevel
                if (L2.eq.0) then
                 write(*,*)dMeanlevel(nv,lv), ' not found for V'
                 stop
                end if


         
                write(*,*)'L1 for U, L2 for V=',L1,L2
 
                do igrid = 1,jf
                  Uapoint=rawdata(igrid,:,ID_U,L1)
                  Vapoint=rawdata(igrid,:,ID_V,L2)
           
                  if(igrid.eq.2328) then
                    kk = 1
                  else
                    kk = 0
                  end if

                  call getvrnce(Uapoint,Vapoint, iens,
     +              UUvrnce,VVvrnce,UVcorr,
     +              Usigma,Vsigma,dMeanlevel(nv,lv),kk)

                  derv_dtra(igrid,lv,1) = UUvrnce
                  derv_dtra(igrid,lv,2) = VVvrnce
                  derv_dtra(igrid,lv,3) = UVcorr
                  derv_dtra(igrid,lv,4) = Usigma
                  derv_dtra(igrid,lv,5) = Vsigma
 
                end do
 
              end do

           end if

        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   subroutine getvrnce: compute mean and spread of one dimension array
c   Author: Binbin Zhou
c   Aug. 3, 2005
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                                                                                    
        subroutine getvrnce (u,v,n,
     +    UUvrnce,VVvrnce, UVcorr,
     +    Usigma,Vsigma,p,kk)

         real u(*), v(*)  
         real UUvrnce,VVvrnce,UVcorr,Usigma,Vsigma
         real Umean, Vmean 
         real x,y
         integer n,kk,p

        
   
         Umean = 0.
         Vmean = 0.

         UUvrnce = 0.
         VVvrnce = 0.
         UVcorr = 0.

         do i=1,n
           Umean = Umean + u(i)
           Vmean = Vmean + v(i)
         end do

         Umean = Umean / n
         Vmean = Vmean / n

         do i = 1, n
           UUvrnce = UUvrnce + (u(i)-Umean)**2
           VVvrnce = VUvrnce + (v(i)-Vmean)**2
         end do
 
           UUvrnce = UUvrnce / (n-1)
           VVvrnce = VVvrnce / (n-1)

         Usigma = sqrt (UUvrnce)
         Vsigma = sqrt (VVvrnce)
 
         x = 0.
         y = 0.
         z = 0.
         do i = 1, n
          x = x + (u(i)-Umean)*(v(i)-Vmean)
          y = y + (u(i)-Umean)**2 
          z = z + (v(i)-Vmean)**2
         end do

          y = y /(n-1)
          z = z /(n-1)

         if (y.gt.0.0.and.z.gt.0.0) then
          UVcorr = x / sqrt(y*z)
          UVcorr = UVcorr / (n-1) 
         else
          UVcorr = 0.
         end if


c         if(kk.eq.1) then
c          write(*,'(a2,21f9.2)') 'u=',(u(i),i=1,21)
c          write(*,'(a2,21f9.2)') 'v=',(v(i),i=1,21)
c          write(*,'(a2,21f9.4)') 'w=',(w(i),i=1,21)
c          write(*,'(a2,21f9.2)') 't=',(t(i),i=1,21)
c          write(*,'(a2,21f9.2)') 'q=',(q(i),i=1,21)
c          write(*,*)UUvrnce,VVvrnce,WWvrnce,UWvrnce,
c     +     VWvrnce,WTvrnce,WQvrnce,UVcorr,Usigma,Vsigma
c         end if

         return
         end

