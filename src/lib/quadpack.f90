module quadpack_m
  implicit none

  private
  public :: QUADPACK_ERROR, qags

  interface
    function quadpack_integrand(x) result(f)
      real, intent(in) :: x
      real             :: f
    end function
  end interface

  character(*), parameter :: QUADPACK_ERROR(6) = [ &
    "maximum number of subdivisions achieved",     &
    "roundoff error detected                ",     &
    "extremely bad integrand behaviour      ",     &
    "algorithm does not converge            ",     &
    "integral is probably divergent         ",     &
    "input is invalid                       "      &
  ]

  external :: dqags

contains

  subroutine qags(f, a, b, atol, rtol, result, aerr, neval, ier, limit, last)
    procedure(quadpack_integrand)  :: f
      !! function to integrate
    real,              intent(in)  :: a
      !! lower integration limit
    real,              intent(in)  :: b
      !! upper integration limit
    real,              intent(in)  :: atol
      !! absolute tolerance
    real,              intent(in)  :: rtol
      !! relative tolerance
    real,              intent(out) :: result
      !! output result of integration
    real,    optional, intent(out) :: aerr
      !! optional: output estimate of the modulus of the absolute error, which should equal or exceed abs(i-result)
    integer, optional, intent(out) :: neval
      !! optional: output number of integrand evaluations
    integer, optional, intent(out) :: ier
      !! optional: output error code
    integer, optional, intent(in)  :: limit
      !! optional: maximum number of subintervals (default: 4096)
    integer, optional, intent(out) :: last
      !! optional, output number of subintervals used

    ! local variables
    integer              :: neval_, ier_, limit_, lenw, last_
    integer, allocatable :: iwork(:)
    real                 :: aerr_
    real,    allocatable :: work(:)

    limit_ = 4096
    if (present(limit)) limit_ = limit
    lenw = limit_ * 4

    allocate (iwork(limit_), work(lenw))

    call dqags(f, a, b, atol, rtol, result, aerr_, neval_, ier_, limit_, lenw, last_, iwork, work)

    if (present(aerr )) aerr  = aerr_
    if (present(neval)) neval = neval_
    if (present(ier  )) ier   = ier_
    if (present(last )) last  = last_
  end subroutine

end module