c        parameter (im = 360, jm= 181)
c        real lat(im,jm), long(im,jm)
c        integer im,jm
c
c        call get_LatLon(im,jm,lat,long)
c        do j = 1, jm
c         write(*,'(20f6.1,a4,20f6.1)') 
c     +   ((long(i,j),lat(i,j)), i=1,10),((long(i,j),lat(i,j)), i=351,im)
c        end do
c
c        stop
c        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine get_LatLon: get !To get global-wide LAT and LON by 
c     given number of grids IM (Longtitude) and JM (Latitude)
c     INPUT: im, jm
c     OUTPUT: lat,long
c
c     Author: Binbin Zhou
c     Jan 12, 2007
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_LatLon(im,jm,lat,long)
        integer, intent (IN) :: im, jm
        REAL,dimension(im,jm), intent(INOUT) :: lat,long

        lat = 0.
        long= 0.

         do i = 1, im
           long (i,1:jm) = 1.0*(i-1)
         end do

         do j = 1, jm
           lat(1:im,j) = 90.0 - 1.0*(j-1)
         end do

         return
         end 

     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   subroutine get_vorticity:  compute vorticity from wind components     
c   UWND and VWND for global model. This code was copyed from GFS
c   INPUT: IM,JM,GDLON,GDLAT,UWND,VWND
c   OUTPUT: ABSV (absolute vorticity
c
c   Author: Binbin Zhou
c   Jan 12, 2007
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
      
      SUBROUTINE get_vorticity(IM,JM,GDLON,GDLAT,UWND,VWND,ABSV)
      
      INTEGER, INTENT(IN) :: IM, JM
      REAL,DIMENSION (IM,JM), INTENT (INOUT) :: GDLON,GDLAT,
     +                                          UWND,VWND,ABSV
      REAL,DIMENSION (IM,JM) :: F

      DTR = 3.14159/180.0
      ERAD = 6.376E6               !earth radius
     
      F0 = 2*6.28 / (24.0*3600.0)

      DO I = 1, IM
       DO J = 1, JM
        F(I,J)=F0 * SIN(GDLAT(I,J)*DTR)
       END DO  
      END DO  
  
c Binbin note: Following code get too large vorticity at both north and south poles

      DO J=1, JM
       DO I=2,IM-1 ! need to do different for i=1 and im
        IF(J == 1 .AND. I <= IM/2)then ! Near North pole
           ABSV(I,J)=((VWND(I+1,J)-VWND(I-1,J))
     &          /((GDLON(I+1,J)-GDLON(I-1,J))*DTR)
     &          +(UWND(I+IM/2,J)*COS(GDLAT(I+IM/2,J)*DTR)+UWND(I,J+1)
     &          *COS(GDLAT(I,J+1)*DTR))/((180.-GDLAT(I,J+1)
     &          -GDLAT(I+IM/2,J))*DTR))/(ERAD*COS(GDLAT(I,J)*DTR))
     &          +F(I,J)
       ELSE IF(J == 1 .AND. I > IM/2)then ! Near North pole
           ABSV(I,J)=((VWND(I+1,J)-VWND(I-1,J))
     &          /((GDLON(I+1,J)-GDLON(I-1,J))*DTR)
     &          +(UWND(I-IM/2,J)*COS(GDLAT(I-IM/2,J)*DTR)+UWND(I,J+1)
     &          *COS(GDLAT(I,J+1)*DTR))/((180.-GDLAT(I,J+1)
     &          -GDLAT(I-IM/2,J))*DTR))/(ERAD*COS(GDLAT(I,J)*DTR))
     &          +F(I,J)                
       ELSE IF(J == JM .AND. I <= IM/2)THEN ! Near South Pole
           ABSV(I,J)=((VWND(I+1,J)-VWND(I-1,J))
     &          /((GDLON(I+1,J)-GDLON(I-1,J))*DTR)
     &          +(UWND(I,J-1)*COS(GDLAT(I,J-1)*DTR)+UWND(I+IM/2,J)*
     &          COS(GDLAT(I+IM/2,J)*DTR))/((180.+GDLAT(I,J-1)
     &          +GDLAT(I+IM/2,J))*DTR))/(ERAD*COS(GDLAT(I,J)*DTR))
     &          +F(I,J)
       ELSE IF(J == JM .AND. I > IM/2)THEN ! Near South Pole
           ABSV(I,J)=((VWND(I+1,J)-VWND(I-1,J))
     &          /((GDLON(I+1,J)-GDLON(I-1,J))*DTR)
     &          +(UWND(I,J-1)*COS(GDLAT(I,J-1)*DTR)+UWND(I-IM/2,J)*
     &          COS(GDLAT(I-IM/2,J)*DTR))/((180.+GDLAT(I,J-1)
     &          +GDLAT(I-IM/2,J))*DTR))/(ERAD*COS(GDLAT(I,J)*DTR))
     &          +F(I,J)
       ELSE
           ABSV(I,J)=((VWND(I+1,J)-VWND(I-1,J))
     &          /((GDLON(I+1,J)-GDLON(I-1,J))*DTR)
     &            -(UWND(I,J-1)*COS(GDLAT(I,J-1)*DTR)-UWND(I,J+1)*
     &          COS(GDLAT(I,J+1)*DTR))/((GDLAT(I,J-1)
     &          -GDLAT(I,J+1))*DTR))/(ERAD*COS(GDLAT(I,J)*DTR))
     &          +F(I,J)
        END IF
       END DO        
      END DO   

      ! compute vorticity at I=1 and I=IM
      DO J=1, JM
        I=1
        IF(J == 1)then ! Near North pole
         ABSV(I,J)=((VWND(2,J)-VWND(IM,J))
     &          /(GDLON(2,J)*DTR)
     &          +(UWND(1+IM/2,J)*COS(GDLAT(1+IM/2,J)*DTR)
     &          +UWND(1,J+1)*COS(GDLAT(1,J+1)*DTR))/
     &          ((180.-GDLAT(1,J+1)-GDLAT(1+IM/2,J))*DTR))
     &          /(ERAD*COS(GDLAT(1,J)*DTR))+F(I,J)
        ELSE IF(J == JM)THEN ! Near South Pole
         ABSV(I,J)=((VWND(2,J)-VWND(IM,J))
     &          /(GDLON(2,J)*DTR)
     &          +(UWND(1,J-1)*COS(GDLAT(1,J-1)*DTR)+UWND(1+IM/2,J)*
     &          COS(GDLAT(1+IM/2,J)*DTR))/((180.+GDLAT(1,J-1)
     &          +GDLAT(1+IM/2,J))*DTR))/(ERAD*COS(GDLAT(1,J)*DTR))
     &          +F(I,J)
        ELSE
          ABSV(I,J)=((VWND(2,J)-VWND(IM,J))
     &      /(GDLON(2,J)*DTR)
     &            -(UWND(1,J-1)*COS(GDLAT(1,J-1)*DTR)-UWND(1,J+1)*
     &          COS(GDLAT(1,J+1)*DTR))/((GDLAT(1,J-1)
     &          -GDLAT(1,J+1))*DTR))/(ERAD*COS(GDLAT(1,J)*DTR))
     &          +F(I,J)
        END IF
        I=IM
        IF(J == 1)then ! Near North pole
         ABSV(I,J)=((VWND(1,J)-VWND(IM-1,J))
     &          /((GDLON(1,J)+360.-GDLON(IM-1,J))*DTR)
     &          +(UWND(IM/2,J)*COS(GDLAT(IM/2,J)*DTR)+UWND(IM,J+1)*
     &          COS(GDLAT(IM,J+1)*DTR))/((180.-GDLAT(IM,J+1)
     &          -GDLAT(IM/2,J))*DTR))/(ERAD*COS(GDLAT(IM,J)*DTR))
        ELSE IF(J == JM)THEN ! Near South Pole
         ABSV(I,J)=((VWND(1,J)-VWND(IM-1,J))
     &          /((GDLON(1,J)+360.-GDLON(IM-1,J))*DTR)
     &          +(UWND(IM,J-1)*COS(GDLAT(IM,J-1)*DTR)+UWND(IM/2,J)*
     &          COS(GDLAT(IM/2,J)*DTR))/((180.+GDLAT(IM,J-1)
     &          +GDLAT(IM/2,J))*DTR))/(ERAD*COS(GDLAT(IM,J)*DTR))
        ELSE
         ABSV(I,J)=((VWND(1,J)-VWND(IM-1,J))
     &      /((GDLON(1,J)+360.-GDLON(IM-1,J))*DTR)
     &            -(UWND(IM,J-1)*COS(GDLAT(IM,J-1)*DTR)-UWND(IM,J+1)*
     &          COS(GDLAT(IM,J+1)*DTR))/((GDLAT(IM,J-1)
     &          -GDLAT(IM,J+1))*DTR))/(ERAD*COS(GDLAT(IM,J)*DTR))
        END IF
      END DO          
 
      RETURN
      END



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine vorticity: compute vorticity 10m or high pressure levels
c     first search for index of U, V from the table -> ID_U, ID_V, if found:
c     then search level index from MeanLevel-> L, if found:
c     then compute the vorticity 
c
c     input: nv,rawdata, im,jm,jf, iens
c     output: derv_mn, derv_sp, derv_pr
c     
c     Author: Binbin Zhou, Jan 15, 2007
c      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine vorticity (nv,rawdata, im,jm,jf, iens,  
     +             derv_mn, derv_sp, derv_pr)

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


        INTEGER, intent(IN) :: nv, jf, iens
        REAL,dimension(jf,iens,numvar,maxmlvl),intent(IN) :: rawdata
        REAL,dimension(jf,nderiv,maxmlvl),intent(INOUT) :: derv_mn
        REAL,dimension(jf,nderiv,maxmlvl),intent(INOUT) :: derv_sp
        REAL,dimension(jf,nderiv,maxplvl,maxtlvl),intent(INOUT) :: 
     +               derv_pr

        REAL,dimension(jf) :: U1d, V1D
        REAL,dimension(jf,iens) :: ABSV
        REAL,dimension(im,jm) :: U2D, V2D, lat, lon, absv2D

        real apoint(iens),Uapoint(iens),Vapoint(iens)

          if(dk6(nv).eq.105) then
           ID_U = index_table(k5,k6,33,105,maxvar)      !search index of direct variable in the table
           ID_V = index_table(k5,k6,34,105,maxvar)      !search index of direct variable in the table
          end if

          if(dk6(nv).eq.100) then
           ID_U = index_table(k5,k6,33,100,maxvar)      !search index of direct variable in the table
           ID_V = index_table(k5,k6,34,100,maxvar)      !search index of direct variable in the table
          end if


          if (ID_U .gt.0 .and. ID_V .gt.0 ) then

             !To get global-wide LAT and LON  
             call get_LatLon(im,jm,lat,lon)

             do lv=1,dMlvl(nv)                       !for all pairs of levels

                L1=index_int_array(MeanLevel(ID_U,:),           !search index of 1st level from  
     +                    dMeanlevel(nv,lv),maxmlvl)            !   ID_U-th variable's MeanLevel 
                if (L1.eq.0) then
                 write(*,*) dMeanlevel(nv,lv), ' not found'
                 stop
                end if

                L2=index_int_array(MeanLevel(ID_V,:),           !search index of 1st level from
     +                    dMeanlevel(nv,lv),maxmlvl)            !   ID_V-th variable's MeanLevel
                if (L2.eq.0) then
                 write(*,*) dMeanlevel(nv,lv), ' not found'
                 stop
                end if

                do k=1,iens

                  U1d = rawdata(:,k,ID_U,L1)
                  V1d = rawdata(:,k,ID_V,L2)
                   
                  do j = 1, jm
                  do i = 1, im
                   ij=(j-1)*im + i
                    U2D(i,j) = U1D(ij)
                    V2D(i,j) = V1D(ij)
                   end do
                  end do

                  call get_vorticity(IM,JM,LON,LAT,U2D,V2D,ABSV2D)

                  do j = 1, jm
                   do i = 1, im
                    ij=(j-1)*im + i
                    ABSV(ij,k) = ABSV2D (i,j)
                   end do                    
                  end do
                 
                 end do                  

                 do igrid = 1,jf
                  apoint=absv(igrid,:)
                  
                  call getmean(apoint,iens,amean,aspread)

                  derv_mn(igrid,nv,lv)=amean
                  derv_sp(igrid,nv,lv)=aspread
                 end do

              end do

              do lv=1,dPlvl(nv)

                L1=index_int_array(MeanLevel(ID_U,:),           !search index of 1st level from
     +                     dProblevel(nv,lv),maxmlvl)           !   ID_U-th variable's MeanLevel
                if (L1.eq.0) STOP 313
                                                                                                                                               
                L2=index_int_array(MeanLevel(ID_V,:),           !search index of 1st level from
     +                     dProblevel(nv,lv),maxmlvl)           !   ID_U-th variable's MeanLevel
                if (L2.eq.0) STOP 314


                do k=1,iens
                                                                                                                                                                                                                                                                  
                  U1d = rawdata(:,k,ID_U,L1)
                  V1d = rawdata(:,k,ID_V,L2)
                                                                                                                                                                                                                                                                  
                  do j = 1, jm
                  do i = 1, im
                   ij=(j-1)*im + i
                    U2D(i,j) = U1D(ij)
                    V2D(i,j) = V1D(ij)
                   end do
                  end do
                                                                                                                                                                                                                                                                  
                  call get_vorticity(IM,JM,LON,LAT,U2D,V2D,ABSV2D)
                                                                                                                                                                                                                                                                  
                  do j = 1, jm
                   do i = 1, im
                    ij=(j-1)*im + i
                    ABSV(ij,k) = ABSV2D (i,j)
                   end do
                  end do
                                                                                                                                                                                                                                                                  
                 end do


                do lt = 1, dTlvl(nv)
                              
                 do igrid = 1,jf

                    apoint=absv(igrid,:)

                    if(trim(dop(nv)).ne.'-') then
                     thr1 = dThrs(nv,lt)
                     thr2 = 0.
               call getprob(apoint,iens,thr1,thr2,dop(nv),aprob)
                     derv_pr(igrid,nv,lv,lt)=aprob
                    else
                     if(lt.lt.dTlvl(nv)) then
                       thr1 = dThrs(nv,lt)
                       thr2 = dThrs(nv,lt+1)
               call getprob(apoint,iens,thr1,thr2,dop(nv),aprob)
                       derv_pr(igrid,nv,lv,lt)=aprob
                     end if
                    end if

                 end do

                end do
              end do

           end if

           return
           end
