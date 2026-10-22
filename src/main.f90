!> Main entry point for MESC optimization.
!>
!> Reads optimization parameters from params_val.nml, calls the appropriate
!> model function via the functor interface, and returns the cost function value.
!
program testvmic
    use, intrinsic :: iso_fortran_env, only: stdout => output_unit
    use mesc_precision_module, only: dp
    use mesc_function_module, only: functn

    implicit none

    real(dp) :: fa               !! computed cost
    integer, parameter :: nx = 6 !! number of optimization parameters

    fa = functn(nx)
    write(stdout,*) "cost12", fa

end program testvmic
