 program pmerge1
 ! merged the optimized parameter values from frc into the parameter files and scein for hwsd_soc optimization
 ! read the parameter numbers from the "params1_frc2_f*.txt', and parameter values from "f${i}/scein_${opt}.txt 
 ! and write out the values of parameters into the parameter files for optimization or initial values 
 implicit none
 real*8,  dimension(14)   :: xparam
 real*8,  dimension(6)    :: param
 integer, dimension(6)    :: nparam
 integer n,nx
 real*8   cost
 character*135 chdata
 
    open(11,file='last_nx_frc.txt')
    open(12,file='last_xp_frc.txt')
    open(13,file='params1_hwsd.txt')
    open(21,file='params1.txt')   ! for optimization
    ! get the parameter numbers
    read(11,*) (nparam(n),n=1,6)    
    read(12,*) cost, (param(n),n=1,6)
    xparam(:) = 1.0
    do n=1,6
       nx=nparam(n)
       xparam(nx) = param(n)
    enddo
    print *, 'xparam= ', xparam(:)
    ! write the optimized parameter values into "params1.txt"
    do n=1,5
       read(13,'(A)') chdata
       write(21,301)  chdata
    enddo
    write(21,302) xparam(1:14)
    write(21,303)
    close(11)
    close(12)
    close(13)
    close(21)    
301 format(a135) 
302 format(14(f6.4,1x))
303 format('  1    10    7     3   14    5  ')  
  end program pmerge1  
     
    
    
    