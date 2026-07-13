!> Fortran 2008 kind parameters used throughout MESC
module precision_module
  use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
  implicit none
  integer, parameter :: r_2 = SELECTED_REAL_KIND(12, 60) !! default real kind (12-digit precision, exponent range 60)
end module precision_module
