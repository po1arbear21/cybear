#include "../util/macro.f90.inc"

module input_src_m

  use error_m,      only: assert_failed
  use math_m,       only: PI
  use bin_search_m, only: bin_search, BS_LESS

  implicit none

  private
  public input_src
  public const_src
  public polygon_src
  public periodic_src
  public sine_src
  public harmonic_src
  public periodic_polygon_src

  type, abstract :: input_src
    !! input source: provide time-dependent values for input variables
  contains
    procedure(input_src_get), deferred :: get
  end type

  abstract interface
    subroutine input_src_get(this, t, y)
      import input_src
      class(input_src), intent(in)  :: this
      real,             intent(in)  :: t
        !! time
      real,             intent(out) :: y(:)
        !! output y(t)
    end subroutine
  end interface

  type, extends(input_src) :: const_src
    !! constant input source
    real, allocatable :: const(:)
      !! output value
  contains
    procedure :: init => const_src_init
    procedure :: get  => const_src_get
  end type

  type, extends(input_src) :: polygon_src
    !! polygonal input source (linear interpolation between predefined nodes)
    real, allocatable :: t(:)
      !! time coordinates of nodes
    real, allocatable :: y(:,:)
      !! y(t) at nodes
  contains
    procedure :: init => polygon_src_init
    procedure :: get  => polygon_src_get
  end type

  type, abstract, extends(input_src) :: periodic_src
    !! periodic inpuit source
    real :: freq
      !! frequency
  contains
    procedure :: periodic_src_init
  end type

  type, extends(periodic_src) :: sine_src
    !! single sine input source
    real, allocatable :: ampl(:)
      !! sine wave amplitude
    real, allocatable :: phase(:)
      !! phase shift
  contains
    procedure :: init => sine_src_init
    procedure :: get  => sine_src_get
  end type

  type, extends(periodic_src) :: harmonic_src
    !! harmonic input source
    real, allocatable :: c(:,:)
      !! cosine coefficients 0:NH
    real, allocatable :: s(:,:)
      !! sine coefficients 1:NH
  contains
    procedure :: init => harmonic_src_init
    procedure :: get  => harmonic_src_get
  end type

  type, extends(periodic_src) :: periodic_polygon_src
    !! periodic polygonal input source
    real, allocatable :: t(:)
      !! time coordinates of nodes
    real, allocatable :: y(:,:)
      !! y(t) at nodes
  contains
    procedure :: init => periodic_polygon_src_init
    procedure :: get  => periodic_polygon_src_get
  end type

contains

  subroutine const_src_init(this, const)
    !! initialize const input source
    class(const_src), intent(out) :: this
    real,             intent(in)  :: const(:)
      !! constant value

    this%const = const
  end subroutine

  subroutine const_src_get(this, t, y)
    !! get value for time t
    class(const_src), intent(in)  :: this
    real,             intent(in)  :: t
      !! time
    real,             intent(out) :: y(:)
      !! return y(t)

    IGNORE(t)
    y = this%const
  end subroutine

  subroutine polygon_src_init(this, t, y)
    !! initialize polygonal input source
    class(polygon_src), intent(out) :: this
    real,               intent(in)  :: t(:)
      !! time coordinates of nodes (must be sorted in ascending order)
    real,               intent(in)  :: y(:,:)
      !! y(t) at nodes; N_val x N_t

    ASSERT(size(y,2) == size(t))

    this%t = t
    this%y = y
  end subroutine

  subroutine polygon_src_get(this, t, y)
    !! get value for time t
    class(polygon_src), intent(in)  :: this
    real,               intent(in)  :: t
      !! time
    real,               intent(out) :: y(:)
      !! return y(t)

    integer :: i
    real    :: w1, w2

    if (t <= this%t(1)) then ! left of left-most node
      y = this%y(:,1)
    elseif (t >= this%t(size(this%t))) then ! right of right-most node
      y = this%y(:,size(this%t))
    else ! somewhere in the middle
      ! get interval i so that t in [t(i):t(i+1)]
      i = bin_search(this%t, t, BS_LESS)

      ! linear weights
      if (this%t(i+1) == this%t(i)) then
        ! avoid division by zero
        w1 = 0.5
        w2 = 0.5
      else
        w2 = (t - this%t(i)) / (this%t(i+1) - this%t(i))
        w1 = 1 - w2
      end if

      ! result
      y = w1 * this%y(:,i) + w2 * this%y(:,i+1)
    end if
  end subroutine

  subroutine periodic_src_init(this, freq)
    !! initialize periodic input source
    class(periodic_src), intent(out) :: this
    real,                intent(in)  :: freq
      !! frequency

    this%freq = freq
  end subroutine

  subroutine sine_src_init(this, freq, ampl, phase)
    !! initialize sine input source
    class(sine_src), intent(out) :: this
    real,            intent(in)  :: freq
      !! frequency
    real,            intent(in)  :: ampl(:)
      !! amplitude
    real, optional,  intent(in)  :: phase(:)
      !! optional phase shift in rad (default: 0)

    ! init base
    call this%periodic_src_init(freq)

    this%ampl  = ampl
    allocate (this%phase(size(ampl)), source = 0.0)
    if (present(phase)) this%phase = phase
  end subroutine

  subroutine sine_src_get(this, t, y)
    !! get value for time t
    class(sine_src), intent(in)  :: this
    real,            intent(in)  :: t
      !! time
    real,            intent(out) :: y(:)
      !! return y(t)

    y = this%ampl * sin(2 * PI * this%freq * t + this%phase)
  end subroutine

  subroutine harmonic_src_init(this, freq, c, s)
    !! initialize harmonic input source
    class(harmonic_src), intent(out) :: this
    real,                intent(in)  :: freq
      !! fundamental frequency
    real,                intent(in)  :: c(:,0:)
      !! cosine coefficients
    real,                intent(in)  :: s(:,:)
      !! sine coefficients

    ASSERT(size(c,1) == size(s,1))
    ASSERT(ubound(c,2) == ubound(s,2))

    ! init base
    call this%periodic_src_init(freq)

    this%c = c
    this%s = s
  end subroutine

  subroutine harmonic_src_get(this, t, y)
    !! get value for time t
    class(harmonic_src), intent(in)  :: this
    real,                intent(in)  :: t
      !! time
    real,                intent(out) :: y(:)
      !! return y(t)

    integer :: i

    y = this%c(:,0)
    do i = 1, ubound(this%c,2)
      y = y + this%c(:,i) * cos(i * 2 * PI * this%freq * t) &
        &   + this%s(:,i) * sin(i * 2 * PI * this%freq * t)
    end do
  end subroutine

  subroutine periodic_polygon_src_init(this, freq, t, y)
    !! initialize polygonal input source
    class(periodic_polygon_src), intent(out) :: this
    real,                        intent(in)  :: freq
      !! frequency
    real,                        intent(in)  :: t(:)
      !! time coordinates of nodes (must be sorted in ascending order, must lie in interval [0; 1/freq])
    real,                        intent(in)  :: y(:,:)
      !! y(t) at nodes; N_val x N_t

    ASSERT(size(y,2) == size(t))
    ASSERT(t(1) >= 0.0)
    ASSERT(t(size(t)) <= 1.0 / freq)

    ! init base
    call this%periodic_src_init(freq)

    this%t = t
    this%y = y
  end subroutine

  subroutine periodic_polygon_src_get(this, t, y)
    !! get value for time t
    class(periodic_polygon_src), intent(in)  :: this
    real,                        intent(in)  :: t
      !! time
    real,                        intent(out) :: y(:)
      !! return y(t)

    integer :: i1, i2, n
    real    :: tperiod, t_, t1, t2, w1, w2

    ! reduce range so that t in [0; 1/freq]
    tperiod = 1.0 / this%freq
    if ((t >= 0) .and. (t <= tperiod)) then
      t_ = t
    else
      t_ = t - floor(t * this%freq) * tperiod
    end if

    ! get interval i1:i2 so that t in [t(i1):t(i2)], extend periodically
    n = size(this%t)
    if (t_ <= this%t(1)) then
      i1 = n
      i2 = 1
      t1 = this%t(n) - tperiod
      t2 = this%t(1)
    elseif (t_ >= this%t(n)) then
      i1 = n
      i2 = 1
      t1 = this%t(n)
      t2 = this%t(1) + tperiod
    else
      i1 = bin_search(this%t, t_, BS_LESS)
      i2 = i1 + 1
      t1 = this%t(i1)
      t2 = this%t(i2)
    end if

    ! linear weights
    if (t2 == t1) then
      ! avoid division by zero
      w1 = 0.5
      w2 = 0.5
    else
      w2 = (t_ - t1) / (t2 - t1)
      w1 = 1 - w2
    end if

    ! result
    y = w1 * this%y(:,i1) + w2 * this%y(:,i2)
  end subroutine

end module
