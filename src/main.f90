program main
  use scheme
  use constant
  implicit none

  call scheme_init()


! contains

!   function uExact (i)
!     real(kind=WP) , dimension(scheme_n) :: uExact, tmp
!     integer, intent(in) :: i
!     real(kind=WP) :: t

!     t = scheme_dt * real( i, WP )
!     if ( scheme_flag == 2 ) then
!        tmp = sin( PI * scheme_grid_size ) / ( PI * scheme_grid_size )
!        uExact = tmp * sin( 2.d0 * PI * ( scheme_grid_node - scheme_a * t ) )
!     else
!        uExact = sin( 2.d0 * PI * (scheme_grid_node - scheme_a * t ))
!     end if

!   end function uExact

!   subroutine getNumericalU ( un )
!     implicit none

!     real(kind=WP), intent(inout), dimension(scheme_n) :: un

!     do i = 1,scheme_n
!        if ( scheme_flag == -2 ) then
!           if ( i-1 < 1 ) then
!              im1 = i - 1 + scheme_n
!           else
!              im1 = i - 1
!           endif
!           if ( i+1 > scheme_n ) then
!              ip1 = i + 1 - scheme_n
!           else
!              ip1 = i + 1
!           end if
!           if ( i+2 > scheme_n ) then
!              ip2 = i + 2 - scheme_n
!           else
!              ip2 = i + 2
!           end if
!           tmp1 = ( scheme_u(i) + scheme_u(ip1) ) * 5.d-1
!           if ( scheme_a > 0 ) then
!              tmp2 = ( scheme_u(im1) - 2.d0 * scheme_u(i) + scheme_u(ip1) ) / 6.d0
!           else
!              tmp2 = ( scheme_u(i) - 2.d0 * scheme_u(ip2) + scheme_u(ip2) ) / 6.d0
!           end if
!           un(i) = tmp1 + tmp2
!        else
!           un(i) = scheme_u(i)
!        end if
!     end do

!   end subroutine getNumericalU

end program main
