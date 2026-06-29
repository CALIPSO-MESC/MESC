program testvmic
    use precision_module, only: dp
    use function_module, only: functn
    implicit none
    real(dp) :: xparam(16),xcost,fa
    integer :: nx

    nx=6
    open(20,file='params_val.txt')
    read(20,*) xcost,xparam(1:nx)
    close(20)

    fa = functn(nx,xparam)
    print *, 'cost12', xcost,xparam(1:nx),fa
    stop
end program testvmic
