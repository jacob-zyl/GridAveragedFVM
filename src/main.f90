program main
  use scheme
  use constant
  implicit none

  integer :: i

  call scheme_init()

  do i = 1, scheme_maxStep
     call scheme_update ()
  end do

end program main
