program main
  use scheme2
  use constant
  implicit none

  integer :: i

  call scheme_init()

  call scheme_writeToFile ( 0 )

  do i = 1, scheme_maxStep
     call scheme_update ()
  end do

  call scheme_writeToFile ( 1 )
  write(*,*) scheme_calculateError( scheme_uExact() ), scheme_dx

end program
