C  raw data
       use grib_mod
       real,allocatable,dimension(:) :: pi, dp1         !jf: total precip  
       real,allocatable,dimension(:) :: cnp1, lgp1      !jf: conv and large scale precip       
       real,allocatable,dimension(:,:) :: t2m 
       real,allocatable,dimension(:) :: dt2m
       real,allocatable,dimension(:,:,:) :: tpr
       real,allocatable,dimension(:,:) :: dtpr

       integer iyr,imon,idy,ihr
       character*50 gdss(400)
       integer GRIBID, kgdss(200), lengds,im,jm,km,jf
       character*40 filename1,filename2,fname
       integer p(10)
       integer iunit,ounit,cyc, hr1,hr2
       type(gribfield) :: gfld
       character*3 fog
       integer fhr
 

         data (p(i),i=1,10)
     +  /1000,975,950,925,900,875,850,825,800,775/


       read(*,*) filename1,filename2,cyc,hr1,hr2,GRIBID,fog,fhr     !hr1: previous fcst hour
       write(*,*)'read ',filename1,filename2,cyc,hr1,hr2,GRIBID,fhr

       write(*,*) filename1,filename2,cyc,hr1,hr2,GRIBID,fhr

       if(GRIBID.eq.255) then   !For NARRE 13km RAP grid#130
         im=1799
         jm=1059
         jf=im*jm
       else
         call makgds(GRIBID, kgdss, gdss, lengds, ier)
         im=kgdss(2)
         jm=kgdss(3)
         jf=kgdss(2)*kgdss(3)
       end if

       write(*,*) 'im jm jf=',im,jm,jf

       allocate(pi(jf))
       allocate(dp1(jf))
       allocate(cnp1(jf))
       allocate(lgp1(jf)) 

       if(fog.eq.'yes') then
        allocate(t2m(jf,2))
        allocate(dt2m(jf))
        allocate(tpr(jf,10,2))
        allocate(dtpr(jf,10))
       end if


        iunit1=10
        iunit2=20

        pi=0.0

        call baopenr(iunit1,filename1,ierr)
         write(*,*) 'open ', trim(filename1), 'ierr=',ierr
         if (ierr.gt.0 ) stop
        call baopenr(iunit2,filename2,ierr)
         write(*,*) 'open ', trim(filename2), 'ierr=',ierr
         if (ierr.gt.0 ) stop
         

ccc  Read Temperature
       if(fog.eq.'yes'.and.fhr.gt.1) then
         jpdtn=0
         jp27=-9999
         call readGB2(iunit1,jpdtn,0,0,103,2,jp27,gfld)
         T2m(:,1)=gfld%fld
         call readGB2(iunit2,jpdtn,0,0,103,2,jp27,gfld)
         T2m(:,2)=gfld%fld

         write(*,*) 'read t2m done'

         do k=1,10
           call readGB2(iunit1,jpdtn,0,0,100,p(k),jp27,gfld)
           Tpr(:,k,1)=gfld%fld
         end do

         do k=1,10
           call readGB2(iunit2,jpdtn,0,0,100,p(k),jp27,gfld)
           Tpr(:,k,2)=gfld%fld
         end do


          write(*,*) 'read tpr done'
        else
          T2m(:,1)=0.
          T2m(:,2)=0.
          Tpr=0.
        end if

ccc  Read APCP
 
        jpdtn=8    !APCP's Product Template# is  4.8 
        jpd1=1
        jpd2=8

        if(filename2(1:3).eq.'rap') then   !RAP members
           jp27=1
    
           !After 2018, RAP drop large-scale APCP, add total APCP
           !call readGB2(iunit2,jpdtn,1,9,1,0,jp27,gfld)
           !write(*,*) 'RAP large APCP getgb2 iret=',iret
           !lgp1(:)=gfld%fld                    !large scale APCP 
           !So directly get total APCP  
            call readGB2(iunit2,jpdtn,1,8,1,0,jp27,gfld)
            write(*,*) 'RAP total APCP getgb2 iret=',iret
            dp1(:)=gfld%fld

            call readGB2(iunit2,jpdtn,1,10,1,0,jp27,gfld)
            write(*,*) 'RAP CONV APCP  getgb2 iret=',iret
            cnp1(:)=gfld%fld                   !Convective APCP
            
            !After 2018, RAP drop large-scale APCP, add total APCP
            !dp1(:)=lgp1(:) + cnp1(:)    !Total APCP 

         else      ! NAMnest members
           if (cyc.eq.0.or.cyc.eq.12) then

               if (hr2.eq.1.or.hr2.eq.13.or.hr2.eq.25 ) then

                 call readGB2(iunit2,jpdtn,1,8,1,0,1,gfld)             
                  dp1(:)=gfld%fld                          !Total APCP
                 call readGB2(iunit2,jpdtn,1,10,1,0,1,gfld)
                  cnp1(:)=gfld%fld                         !Convective APCP

               else
        
                 !need previous fcst hour data
                 call readGB2(iunit1,jpdtn,1,8,1,0,-9999,gfld)
                  dp1(:)=gfld%fld
                 call readGB2(iunit1,jpdtn,1,10,1,0,-9999,gfld)
                  cnp1(:)=gfld%fld

                 jp27=-9999
                 if(hr2.eq.3) then jp27=3
                 if(hr2.eq.6) then jp27=6
                 if(hr2.eq.9) then jp27=9
                 if(hr2.eq.12) then jp27=12
                 if(hr2.eq.15) then jp27=3
                 if(hr2.eq.18) then jp27=6
                 if(hr2.eq.21) then jp27=9
                 if(hr2.eq.24) then jp27=12
                 if(hr2.eq.27) then jp27=3
                 if(hr2.eq.30) then jp27=6
                 if(hr2.eq.33) then jp27=9
                 if(hr2.eq.36) then jp27=12
                 call readGB2(iunit2,jpdtn,1,8,1,0,jp27,gfld)
                  dp1(:)=gfld%fld-dp1(:)
                 call readGB2(iunit2,jpdtn,1,10,1,0,jp27,gfld)
                  cnp1(:)=gfld%fld-cnp1(:)

               end if          
           else      !cyc=06, 18Z

               if( hr2.eq.1.or.hr2.eq.4.or.hr2.eq.7.or.
     +             hr2.eq.10.or.hr2.eq.13.or.hr2.eq.16.or.
     +             hr2.eq.19.or.hr2.eq.22.or.hr2.eq.25.or.
     +             hr2.eq.28.or.hr2.eq.31.or.hr2.eq.34) then
                
                   call readGB2(iunit2,jpdtn,1,8,1,0,1,gfld)
                     dp1(:)=gfld%fld
                   call readGB2(iunit2,jpdtn,1,10,1,0,1,gfld)
                     cnp1(:)=gfld%fld
             
                else

                   call readGB2(iunit1,jpdtn,1,8,1,0,-9999,gfld) !need previous fcst hour data
                      dp1(:)=gfld%fld
                   call readGB2(iunit1,jpdtn,1,10,1,0,-9999,gfld) !need previous fcst hour data
                      cnp1(:)=gfld%fld
                   call readGB2(iunit2,jpdtn,1,8,1,0,-9999,gfld)
                      dp1(:)=gfld%fld-dp1(:)
                   call readGB2(iunit2,jpdtn,1,10,1,0,-9999,gfld)
                      cnp1(:)=gfld%fld-cnp1(:)

                 end if

           end if
         end if 

        call baclose(iunit1,ierr)
        call baclose(iunit2,ierr)
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Compute temperature change hourly or 3 hourly
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if(fog.eq.'yes') then
         dt2m(:)=t2m(:,2)-t2m(:,1)
         do k=1,10
          dtpr(:,k)=tpr(:,k,2)-tpr(:,k,1)
         end do
       end if 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc    Compute 3hr, 6hr accumulated precip
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       !do i=8500,8600
       ! write(*,'(i10,a10,f10.3)') i, ' cnp1,dp1=',cnp1(i),dp1(i)
       ! if(fog.eq.'yes') then
       !   write(*,'(a7,3f10.2)')  '   T2m=',t2m(i,1),t2m(i,2),dt2m(i)
       !   write(*,'(a7,10f10.2)') 'tpr(1)=',(tpr(i,k,1),k=1,10)
       !   write(*,'(a7,10f10.2)') 'tpr(2)=',(tpr(i,k,2),k=1,10)
       !   write(*,'(a7,10f10.2)') '  dtpr=',(dtpr(i,k),k=1,10)
       ! end if
       !end do



cccccc  Then call putgb2 to store the calculated data into a grib2 file
c
c      data structure gfld is re-used for pack data since all are same
c      only gfld%fld and gfld%ipdtmpl(27) are different
c
        ounit=30
        fname=trim(filename2)//'.prcp'
        call baopen(ounit,fname,ierr)
        
        gfld%ipdtnum=8

        !Store APCP
        gfld%ipdtmpl(1)=1
        gfld%ipdtmpl(2)=8     !Total APCP
        gfld%ipdtmpl(10)=1
        gfld%ipdtmpl(12)=0

        gfld%fld(:)=dp1(:) 
        gfld%ipdtmpl(9)=hr2-1
        if(gfld%ipdtmpl(9).lt.0) gfld%ipdtmpl(9)=0
        gfld%ipdtmpl(27)=1
        call putgb2(ounit,gfld,ierr)
        if(ierr.ne.0) write(*,*) 'packing dp1 error'

        gfld%ipdtmpl(2)=10    !Convective APCP
        gfld%fld(:)=cnp1(:)
        call putgb2(ounit,gfld,ierr)
        if(ierr.ne.0) write(*,*) 'packing cnp1 error'

        !Store Temperature variation
        if(fog.eq.'yes') then 
          gfld%fld(:)=dt2m(:)
          gfld%ipdtmpl(1)=0
          gfld%ipdtmpl(2)=0
          gfld%ipdtmpl(10)=103
          gfld%ipdtmpl(12)=2
          gfld%ipdtmpl(9)=hr2-1
          gfld%ipdtmpl(27)=1
          call putgb2(ounit,gfld,ierr)
          if(ierr.ne.0) write(*,*) 'packing dt2m error'

          do k=1,10
            gfld%fld(:)=dtpr(:,k)
            gfld%ipdtmpl(1)=0
            gfld%ipdtmpl(2)=0
            gfld%ipdtmpl(10)=100
            gfld%ipdtmpl(12)=p(k)*100
            gfld%ipdtmpl(9)=hr2-1
            gfld%ipdtmpl(27)=1
            call putgb2(ounit,gfld,ierr)
            if(ierr.ne.0) write(*,*) 'packing dtpr error'
          end do  
         end if

         call baclose(ounit,ierr)   
         write(*,*) 'close ', fname, 'ierr=',ierr

      stop
      end


      subroutine readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,gfld)

        use grib_mod

        type(gribfield) :: gfld 
 
        integer jids(200), jpdt(200), jgdt(200)
        integer jpd1,jpd2,jpdtn
        logical :: unpck=.true. 

        jids=-9999  !array define center, master/local table, year,month,day, hour, etc, -9999 wildcard to accept any
        jpdt=-9999  !array define Product, to be determined
        jgdt=-9999  !array define Grid , -9999 wildcard to accept any

        jdisc=-1    !discipline#  -1 wildcard 
        jgdtn=-1    !grid template number,    -1 wildcard 
        jskp=0      !Number of fields to be skip, 0 search from beginning
        ifile=0

        jpdt(1)=jpd1   !Category #     
        jpdt(2)=jpd2   !Product # under this category     
        jpdt(10)=jpd10
        if(jpd10.eq.100) then
           jpdt(12)=jpd12*100   !pressure level     
        else
           jpdt(12)=jpd12
        end if

        jpdt(27)=jpd27  !Time range (1 hour, 3 hr etc)
        

         call getgb2(igrb2,ifile,jskp,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     +        unpck, jskp1, gfld,iret)

        return
        end 
         
