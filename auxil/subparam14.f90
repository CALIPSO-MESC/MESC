    program substitue
    ! produce a line of 14 parameter values for replacing default values in "params1.txt"
    implicit none
    integer i,nx(6)
    real cost,x(14),xp(6)

      open(1,file='last2.txt')
      open(2,file='p14.txt')
      x(:) = 1.0
      read(1,*) (nx(i),i=1,6)
      read(1,*) cost,(xp(i),i=1,6)
      do i=1,6
         x(nx(i)) =xp(i)
      enddo
      x(12)=1.0
      write(2,21) x(1:14)

      close(1)
      close(2)
21   format(14(f8.4,1x))
    end program substitue
