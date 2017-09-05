program main
  use scheme
  use constant
  implicit none

  integer :: i

  call scheme_init()

  call scheme_writeToFile ( 0 )

  do i = 1, scheme_maxStep
     call scheme_update ()
  end do

  call scheme_writeToFile ( 1 )
  write(*,*) scheme_calculateErrorInfinity( scheme_uExact() ), &
       &     scheme_calculateErrorL1( scheme_uExact() ), scheme_dx
end program
