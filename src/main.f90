!> Main entry point for MESC optimization.
!>
!> Reads optimization parameters from params_val.nml, calls the appropriate
!> model function via the functor interface, and returns the cost function value.
!
program testvmic

    use precision_module, only: dp
    use function_module, only: functn

    implicit none

    real(dp) :: xparam(16)       !! optimization params
    real(dp) :: fa               !! computed cost
    integer, parameter :: nx = 6 !! number of optimization parameters
    namelist /params/ xparam     !! define the namelist

    ! Read from the namelist
    open(20, file="params_val.nml", status="old")
    read(20, nml=params)
    close(20)

    fa = functn(nx, xparam)
    print *, "cost12", xparam(1:nx), fa

end program testvmic
