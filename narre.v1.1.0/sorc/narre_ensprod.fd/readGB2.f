      subroutine readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +     gfld,iret)

        use grib_mod

        type(gribfield),intent(IN) :: gfld
        integer jids(200), jpdt(200), jgdt(200)
        integer igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27
        logical :: unpack=.true.

c        write(*,*) 'igrb2=',igrb2, jpdtn,jpd1,jpd2,jpd10,jpd12

        jids=-9999  !array define center, master/local table, year,month,day, hour, etc, -9999 wildcard to accept any
        jpdt=-9999  !array define Product, to be determined
        jgdt=-9999  !array define Grid , -9999 wildcard to accept any

        jdisc=-1    !discipline#  -1 wildcard 
        jgdtn=-1    !grid template number,    -1 wildcard 
        jskp=0      !Number of fields to be skip, 0 search from beginning

        jpdt(1)=jpd1   !Category #     
        jpdt(2)=jpd2   !Product # under this category     
        jpdt(10)=jpd10 !Product vertical ID      
        jpdt(27)=jpd27

        if(jpd10.eq.100) then
           jpdt(12)=jpd12*100   !pressure level     
        else 
           jpdt(12)=jpd12
        end if

        call getgb2(igrb2,0,jskp,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     +     unpack, jskp1, gfld,iret)

c        if(iret.ne.0) then
c         write(*,*) 'getgb2 error:',iret,' in read ',igrb2,
c     +   ' for var ', jpd1,jpd2,jpd10,jpd12
c        end if

        return
        end
