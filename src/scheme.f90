module scheme
  use constant
  implicit none
  private

  real(kind=WP), parameter :: a = 1.0 ! coefficient in equation
  real(kind=WP), parameter :: L = 1.0 ! length of spatial region to solve
  real(kind=WP), parameter :: T = 1.5 ! length of temporal region to solve

  real(kind=WP), public :: scheme_CFL, scheme_dt, scheme_dx
  integer, public :: scheme_maxStep, scheme_numOfGrid

  ! the primary data of region to solve
  real(kind=WP), public, allocatable, dimension(:) :: scheme_u

  real(kind=WP), public, allocatable, dimension(:) :: scheme_gridSize
  real(kind=WP), public, allocatable, dimension(:) :: scheme_gridNode
  real(kind=WP), public, allocatable, dimension(:) :: scheme_flux

  public :: scheme_init, scheme_update, scheme_calculateError,&
       & scheme_writeToFile, scheme_uExact


contains

  subroutine scheme_init ()
    implicit none
    integer :: i
    real(kind=WP) :: tmp

    !! First, read and set input data
    !====================================
    ! Format of input data file (data.in)
    !====================================
    !0.99 \tab ! CFL
    !10   \tab ! numOfGrid
    !====================================
    open(unit=11, file='data.in', action='read', status='old')
    read(11, *) scheme_CFL
    read(11, *) scheme_numOfGrid
    scheme_dx = L / scheme_numOfGrid
    scheme_dt = scheme_dx * scheme_CFL / a
    scheme_maxStep = floor( T / scheme_dt ) + 1
    scheme_dt = T / scheme_maxStep

    !! Second, allocate arrays
    associate( n => scheme_numOfGrid )
      allocate( scheme_gridSize(-1:n+2), scheme_gridNode(-1:n+2), &
           & scheme_u(-1:n+2), scheme_flux(0:n) )
    end associate

    !! Third, set data of arrays
    scheme_gridSize = scheme_dx
    forall ( i = -1:scheme_numOfGrid+2 )
       scheme_gridNode(i) = ( i - 0.5 ) * scheme_dx
    end forall

    do i = -1,scheme_numOfGrid+2
       tmp = sin( PI * scheme_gridSize(i) ) / ( PI * scheme_gridSize(i) )
       scheme_u(i) = tmp * sin( 2.0 * PI * scheme_gridNode(i) )
    end do
  end subroutine scheme_init

  subroutine scheme_update ()
    implicit none
    ! Update scheme for *ONE* step
    integer :: i
    real(kind=WP), allocatable, dimension(:) :: k1, k2, k3, k4
    real(kind=WP), allocatable, dimension(:) :: u1

    associate(n=>scheme_numOfGrid)
      allocate( u1(1:n) )
      allocate( k1(1:n), k2(1:n), k3(1:n), k4(1:n) )

      ! RK4 method

      u1 = scheme_u(1:n)

      call flux_update()
      forall ( i = 1:n )
         k1(i) = scheme_dt * ( scheme_flux(i-1) - scheme_flux(i) ) / scheme_dx
      end forall

      scheme_u(1:n) = u1 + 0.5 * k1
      call flux_update()

      forall ( i = 1:n )
         k2(i) = scheme_dt * ( scheme_flux(i-1) - scheme_flux(i) ) / scheme_dx
      end forall

      scheme_u(1:n) = u1 + 0.5 * k2
      call flux_update()

      forall ( i = 1:n )
         k3(i) = scheme_dt * ( scheme_flux(i-1) - scheme_flux(i) ) / scheme_dx
      end forall

      scheme_u(1:n) = u1 + k3
      call flux_update()

      forall ( i = 1:n )
         k4(i) = scheme_dt * ( scheme_flux(i-1) - scheme_flux(i) ) / scheme_dx
      end forall

      scheme_u(1:n) = u1 + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0

      ! End of RK4

    end associate
    call scheme_boudaryCondition ()

  contains

    function kappa (a)
      real(kind=WP) :: kappa
      real(kind=WP), intent(in) :: a

      if ( abs(a) > epsilon(a) ) then
         kappa = a / abs(a)
      else
         kappa = a / (a*a + epsilon(a)*epsilon(a)) / (2 * epsilon(a))
      end if
    end function kappa

    subroutine flux_update ()
      implicit none
      integer :: i
      real(kind=WP) :: tmp1, tmp2, tmp3

      call scheme_boudaryCondition ()
      do i = 1, scheme_numOfGrid
         tmp1 = scheme_u(i) + scheme_u(i+1)
         tmp2 = scheme_u(i-1) + scheme_u(i+2)
         tmp3 = scheme_u(i+2) - 3.0*scheme_u(i+1) + 3.0*scheme_u(i) - scheme_u(i-1)
         scheme_flux(i) = 7.0*a/12.0*tmp1 - tmp2*a/12.0 + kappa(a)*a/12.0*tmp3
      end do
      scheme_flux(0) = scheme_flux(scheme_numOfGrid)
    end subroutine flux_update

  end subroutine scheme_update

  function scheme_calculateError (ue) result (error)
    implicit none
    real(kind=WP), intent(in), dimension(scheme_numOfGrid) :: ue
    real(kind=WP) :: error

    associate( n => scheme_numOfGrid )
      error = maxval( abs( ue(1:n) - scheme_u(1:n) ) )
    end associate
  end function scheme_calculateError

  subroutine scheme_boudaryCondition ()
    implicit none
    ! This is periodic boundary condition

    associate(n=>scheme_numOfGrid)
      scheme_u(0) = scheme_u(n)
      scheme_u(-1) = scheme_u(n-1)
      scheme_u(n+1) = scheme_u(1)
      scheme_u(n+2) = scheme_u(2)
    end associate

  end subroutine scheme_boudaryCondition

  subroutine scheme_writeToFile ( n )
    implicit none

    integer, intent(in) :: n
    integer :: i
    character( len = 2 ) :: cTemp

    write(cTemp, '(i2)') n

    open(unit=101, file='data' // trim(adjustl(cTemp)) // '.out',&
         & action='write', status='new')
    write(101, *) (scheme_u(i), i = 1, scheme_numOfGrid)
    close(101)
  end subroutine scheme_writeToFile

  function scheme_uExact () result ( u )
    real(kind=WP), dimension(scheme_numOfGrid) :: u, tmp

    associate( n => scheme_numOfGrid )
      tmp = sin( PI * scheme_gridSize(1:n) ) / ( PI * scheme_gridSize(1:n) )
      u = tmp * sin( 2.0 * PI * ( scheme_gridNode(1:n) - a * T ) )
    end associate
  end function scheme_uExact

end module scheme
