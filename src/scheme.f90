module scheme
  use constant
  implicit none
  private

  real(kind=WP), parameter :: a = 1.0 ! coefficient in equation
  real(kind=WP), parameter :: L = 3.0 ! length of spatial region to solve
  real(kind=WP), parameter :: T = 5.0 ! length of temporal region to solve

  real(kind=WP), public :: scheme_CFL, scheme_dt, scheme_dx
  integer, public :: scheme_maxStep, scheme_numOfGrid

  ! the primary data of region to solve
  real(kind=WP), public, allocatable, dimension(:) :: scheme_u

  real(kind=WP), public, allocatable, dimension(:) :: scheme_gridSize
  real(kind=WP), public, allocatable, dimension(:) :: scheme_gridNode
  real(kind=WP), public, allocatable, dimension(:) :: scheme_flux

  public :: scheme_init, scheme_update, scheme_calculateError


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
    allocate( scheme_gridSize(-1:scheme_numOfGrid+2), &
         & scheme_gridNode(-1:scheme_numOfGrid+2), &
         & scheme_u(-1:scheme_numOfGrid+2), &
         & scheme_flux(0:scheme_numOfGrid) )

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
    !! something

  end subroutine scheme_update

  subroutine scheme_calculateError (arg)
    real, intent(in) :: arg

    
  end subroutine scheme_calculateError





end module scheme
