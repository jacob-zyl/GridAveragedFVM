! Created by Jacob Zeyu LIU (liuzeyu271828@gmail.com)
!
! Must be compiled with -fdefault-real-8 with gfortran or
! -r8 in ifort.

module constant
  implicit none

  integer, parameter :: SGL = selected_real_kind(p=6)
  integer, parameter :: DBL = selected_real_kind(p=13)
  integer, parameter :: WP = DBL
  real(kind=WP), parameter :: E = exp(1.0), PI = 4.0*atan(1.0)

end module constant

