module random_m
  use iso_c_binding
  implicit none

  private
  public :: random

  type, bind(c) :: pcg64_state
    integer(kind=c_int64_t) :: state(2)
    integer(kind=c_int64_t) :: inc(2)
  end type

  type random
    !! random number generator interface (pcg64)
    type(pcg64_state) :: state
      !! internal state
  contains
    procedure :: init       => random_init_
    procedure :: next_int   => random_next_int
    procedure :: next_ints  => random_next_ints
    procedure :: next_real  => random_next_real
    procedure :: next_reals => random_next_reals
  end type

  ! interfaces to pcg64 routines
  interface
    subroutine pcg64_srandom(rng, initstate_h, initstate_l, initseq_h, initseq_l) bind(c)
      import :: pcg64_state, c_int64_t
      type(pcg64_state),              intent(out) :: rng
      integer(kind=c_int64_t), value, intent(in)  :: initstate_h
      integer(kind=c_int64_t), value, intent(in)  :: initstate_l
      integer(kind=c_int64_t), value, intent(in)  :: initseq_h
      integer(kind=c_int64_t), value, intent(in)  :: initseq_l
    end subroutine

    function pcg64_random(rng) result(x) bind(c)
      import :: pcg64_state, c_int64_t
      type(pcg64_state), intent(inout) :: rng
      integer(kind=c_int64_t)          :: x
    end function

    subroutine pcg64_random_n(rng, n, x) bind(c)
      import :: pcg64_state, c_int64_t
      type(pcg64_state),              intent(inout) :: rng
      integer(kind=c_int64_t), value, intent(in)    :: n
      integer(kind=c_int64_t),        intent(out)   :: x(*)
    end subroutine

    function pcg64_boundedrand(rng, bound) result(x) bind(c)
      import :: pcg64_state, c_int64_t
      type(pcg64_state),              intent(inout) :: rng
      integer(kind=c_int64_t), value, intent(in)    :: bound
      integer(kind=c_int64_t)                       :: x
    end function

    subroutine pcg64_boundedrand_n(rng, bound, n, x) bind(c)
      import :: pcg64_state, c_int64_t
      type(pcg64_state),              intent(inout) :: rng
      integer(kind=c_int64_t), value, intent(in)    :: bound
      integer(kind=c_int64_t), value, intent(in)    :: n
      integer(kind=c_int64_t),        intent(out)   :: x(*)
    end subroutine

    function pcg64_random_real(rng) result(x) bind(c)
      import :: pcg64_state, c_double
      type(pcg64_state), intent(inout) :: rng
      real(kind=c_double)              :: x
    end function

    subroutine pcg64_random_real_n(rng, n, x) bind(c)
      import :: pcg64_state, c_int64_t, c_double
      type(pcg64_state),              intent(inout) :: rng
      integer(kind=c_int64_t), value, intent(in)    :: n
      real(kind=c_double),            intent(out)   :: x(*)
    end subroutine

    subroutine pcg64_advance(rng, delta_h, delta_l) bind(c)
      import :: pcg64_state, c_int64_t
      type(pcg64_state),              intent(inout) :: rng
      integer(kind=c_int64_t), value, intent(in)    :: delta_h
      integer(kind=c_int64_t), value, intent(in)    :: delta_l
    end subroutine
  end interface

contains

  subroutine random_init_(this, seed, seq)
    !! initialize random number generator (applies seed and selects sequence)
    class(random),   intent(out) :: this
    integer(kind=8), intent(in)  :: seed(2)
      !! 128-bit seed
    integer(kind=8), intent(in)  :: seq(2)
      !! 128-bit sequence selector (give each thread its own sequence)

    call pcg64_srandom(this%state, seed(1), seed(2), seq(1), seq(2))
  end subroutine

  function random_next_int(this, bound) result(x)
    !! get next equidistributed integer with 0 <= x < bound
    class(random),             intent(inout) :: this
    integer(kind=8), optional, intent(in)    :: bound
      !! optional bound
    integer(kind=8)                          :: x
      !! return random number

    if (present(bound)) then
      x = pcg64_boundedrand(this%state, bound)
    else
      x = pcg64_random(this%state)
    end if
  end function

  subroutine random_next_ints(this, x, bound)
    !! get next equidistributed integers with 0 <= x < bound
    class(random),             intent(inout) :: this
    integer(kind=8),           intent(out)   :: x(:)
      !! output random numbers
    integer(kind=8), optional, intent(in)    :: bound
      !! optional bound

    if (present(bound)) then
      call pcg64_boundedrand_n(this%state, bound, size(x, kind=8), x)
    else
      call pcg64_random_n(this%state, size(x, kind=8), x)
    end if
  end subroutine

  function random_next_real(this) result(x)
    !! get next equidistributed real with 0.0 <= x <= 1.0
    class(random),     intent(inout) :: this
    real                             :: x
      !! return random number

    x = pcg64_random_real(this%state)
  end function

  subroutine random_next_reals(this, x)
    !! get next equidistributed reals with 0 <= x <= 1
    class(random), intent(inout) :: this
    real,          intent(out)   :: x(:)
      !! output random numbers

    call pcg64_random_real_n(this%state, size(x, kind=8), x)
  end subroutine

end module