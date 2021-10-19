C  raw data
       use grib_mod
       real,allocatable,dimension(:,:) :: dp1,dp3,dp6 !jf,12: total precip  

       integer iyr,imon,idy,ihr
       character*50 gdss(400)
       integer IENS, kgdss(200), lengds,im,jm,km,jf
       character*40 filename(100),fname,output
       integer dt(100), ifhr(100),p(14),IGRID
       integer iunit,ounit
       type(gribfield) :: gfld
       character*2 fhr(100),cyc 
       character*3 GRID,mdl

       data (fhr(i),i=1,12)
     + /'01','02','03','04','05','06','07','08','09','10',
     +  '11','12'/

       data (ifhr(i),i=1,12)
     + /1,2,3,4,5,6,7,8,9,10,11,12/

cc     RAP has one-hour accumu precip, so only one file is used
cc     NAM has no one-hour accumu precip, so two files are needed
       read(*,*) IGRID, mdl,cyc 

       if(IGRID.eq.255) then   !For NARRE 13km RAP grid#130
         im=1799
         jm=1059
         jf=im*jm
       else
         call makgds(IGRID, kgdss, gdss, lengds, ier)
         im=kgdss(2)
         jm=kgdss(3)
         jf=kgdss(2)*kgdss(3)
       end if

       if(IGRID.eq.130) then
         GRID='130'
       else if (IGRID.eq.242) then
        GRID='242'
       end if

       write(*,*) 'jf=',jf

       allocate(dp1(jf,12))
       allocate(dp3(jf,12))
       allocate(dp6(jf,12))

       do 1000 i=1,12    !fhr(i):00 03 06 09, ....or 00 01 02 03 04 .... 

        filename(i)='prcip.t'//cyc//'z.'//mdl//'.f'//fhr(i)

        iunit=10+i
        call baopenr(iunit, filename(i),ierr)
        write(*,*) 'open ', filename(i) , 'ierr=',ierr
        if (ierr.ne.0) then
          dp1(:,i)=0.0
          call baclose(iunit,ierr)
          goto 1000
        end if

ccc  Read APCP
 
        jpdtn=8    !APCP's Product Template# is  4.8 
        jpd1=1
        jpd2=8
        jp27=1
        call readGB2(iunit,jpdtn,1,8,1,0,jp27,gfld)
         numpts=gfld%ndpts
         write(*,*) ' ndpts read =',  numpts
         if (numpts.eq.0) then
          dp1(:,i)=0.0
          call baclose(iunit,ierr)
          write(*,*) 'close ', filename(i) , 'ierr=',ierr
         else
           dp1(:,i)=gfld%fld                    
           call baclose(iunit,ierr)
           write(*,*) 'close ', filename(i) , 'ierr=',ierr
         end if
1000  continue

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc    Compute 3hr, 6hr, 12hr and 24hr accumulated precip
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       dp3=0.0
       dp6=0.0

        do i=1,2      
         dp3(:,i)=-1.0
        end do
        do i=3,12    
         dp3(:,i)=dp1(:,i)+dp1(:,i-1)+dp1(:,i-2)
        end do  

        do i=1,5        
         dp6(:,i)=-1.0
        end do
        do i=6,12 
         dp6(:,i)=dp1(:,i)+dp1(:,i-1)+dp1(:,i-2)+dp1(:,i-3)+
     +            dp1(:,i-4)+dp1(:,i-5)
        end do

      

cccccc  Then call putgb2 to store the calculated data into a grib2 file
c
c      data structure gfld is re-used for pack data since all are same
c      only gfld%fld and gfld%ipdtmpl(27) are different
c
       do 2000 i=1,12     !don't consider f00 

       output='hr3_'//filename(i)

        ounit=100+12
        call baopen(ounit,output,ierr)

        gfld%fld(:)=dp3(:,i)
        gfld%ipdtmpl(9)=ifhr(i)-3       
        if(gfld%ipdtmpl(9).lt.0) gfld%ipdtmpl(9)=0
        gfld%ipdtmpl(27)=3       
        call putgb2(ounit,gfld,ierr)

        gfld%fld(:)=dp6(:,i)
        gfld%ipdtmpl(9)=ifhr(i)-6
        if(gfld%ipdtmpl(9).lt.0) gfld%ipdtmpl(9)=0
        gfld%ipdtmpl(27)=6       
        call putgb2(ounit,gfld,ierr)

        call baclose(ounit,ierr)   
          write(*,*) 'close ', output, 'ierr=',ierr

2000   continue

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
         
