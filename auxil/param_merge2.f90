 program pmerge2
 ! merged the optimized parameter values from hwsd_soc into the parameter files and scein for global optimization
 ! read the parameter numbers from the "params1_hwsd2_${case}.txt', and parameter values from "${case}/scein_${opt}.txt
 ! and write out the values of parameters into the parameter files for optimization or initial values
 implicit none
 real(8),  dimension(14)   :: xparam
 real(8),  dimension(6)    :: param
 real(8),  dimension(3)    :: y
 integer, dimension(6)    :: nparam
 integer :: n,nx,ny
 real(8)   :: cost
 character(len=135) :: chdata
 real :: y0,ymin,ymax,y1,y2,y3
 character :: data*30, data1*30,data2*20
    ! input files
    open(11,file='last2_hwsd_param.txt')       ! last two lines in the hwsd_soc parameter file
    open(12,file='last1_hwsd_sceout.txt')      ! last line of sceout from hwsd_soc optimization
    open(13,file='params1_global.txt')         ! the global parameter file to be updated
    open(14,file='last1_global_sceout.txt')    ! last estiamets of global parameter optimization
    open(2,file='scein.dat')                   ! optimization setting

    ! output files
    open(21,file='params1.txt')                ! for optimization setting to be updated
    open(3,file='scein_copy.dat')              ! optimization setting with updated priors

    ! get the parameter numbers
    read(11,*) (xparam(n),n=1,14)
    read(11,*) (nparam(n),n=1,6)
    read(12,*) cost, (param(n),n=1,6)
    xparam(:) = 1.0
    do n=1,6
       nx=nparam(n)
       xparam(nx) = param(n)
    end do
    print *, 'xparam= ', xparam(:)
    ! write the optimized parameter values into "params1.txt"
    do n=1,9
       read(13,'(A)') chdata
       write(21,301)  chdata
    end do
    write(21,302) xparam(1:14)
    write(21,303)

    ! update the initial values in scein.dat
    ny=3
    y(1) = xparam(1)
    y(2) = xparam(10)
    y(3) = xparam(7)

    cost=0.0;y1=0.0;y2=0.0;y(3)=0.0
    read(14,*,end=15) cost,y1,y2,y3
15  if(cost>0.0) then
      y(1) =y1
      y(2) =y2
      y(3) =y3
    end if

    do n=1,2
       read(2,21) data
       write(3,21) data
    end do
    do n=1,ny
       read(2,22) data1,data2
       read(data1,*) y0,ymin,ymax
       write(3,31) y(n),ymin,ymax,data2
    end do

    close(11)
    close(12)
    close(13)
    close(14)
    close(21)
    close(2)
    close(3)

301 format(a135)
302 format(14(f8.4,1x))
303 format('  1    10    7    ')

21  format(a30)
22  format(a30,a20)
31  format(3(f8.2,2x),a20)
  end program pmerge2




