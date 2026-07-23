    program substitue
    ! use the last unfinished SCEUA parameter values as the initial values for next SCEUA
    implicit none
    integer :: nx,ny,n
    real :: y0,ymin,ymax
    character :: data*30, data1*30,data2*20


    real, dimension(:), allocatable  :: varx,vary,y

      nx=1; ny=3
      allocate(varx(nx),vary(ny),y(ny))

      open(1,file='last7.txt')
      open(2,file='scein.dat')
      open(3,file='scein_copy.dat')
      read(1,*) varx(1:nx),y(1:ny)
      do n=1,2
         read(2,201) data
         write(3,201) data
      end do
      do n=1,ny
         read(2,202) data1,data2
         read(data1,*) y0,ymin,ymax
         write(3,301) y(n),ymin,ymax,data2
      end do
      deallocate(varx,vary)
      close(1)
      close(2)
      close(3)
201   format(a30)
202   format(a30,a20)
301   format(3(f8.2,2x),a20)

    end program substitue
