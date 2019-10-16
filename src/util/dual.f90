module dual_m
  use high_precision_m
  implicit none

  type dual
    !! dual number type for automatic differentiation (forward mode)

    integer           :: n
      !! number of base variables
    real              :: x
      !! value of expression
    real, allocatable :: dx(:)
      !! derivatives wrt. base variables
  contains
    procedure :: init => dual_init
  end type

  interface operator(+)
    module procedure :: dual_add_dual
    module procedure :: dual_add_real
    module procedure :: real_add_dual
  end interface

  interface operator(-)
    module procedure :: dual_neg
    module procedure :: dual_sub_dual
    module procedure :: dual_sub_real
    module procedure :: real_sub_dual
  end interface

  interface operator(*)
    module procedure :: dual_mul_dual
    module procedure :: dual_mul_real
    module procedure :: real_mul_dual
  end interface

  interface operator(/)
    module procedure :: dual_div_dual
    module procedure :: dual_div_real
    module procedure :: real_div_dual
  end interface

  interface operator(**)
    module procedure :: dual_pow_dual
    module procedure :: dual_pow_real
    module procedure :: dual_pow_int
    module procedure :: real_pow_dual
  end interface

  intrinsic abs
  interface abs
    module procedure :: dual_abs
  end interface

  intrinsic sqrt
  interface sqrt
    module procedure :: dual_sqrt
  end interface

  intrinsic exp
  interface exp
    module procedure :: dual_exp
  end interface

  intrinsic log
  interface log
    module procedure :: dual_log
  end interface

  intrinsic sin
  interface sin
    module procedure :: dual_sin
  end interface

  intrinsic cos
  interface cos
    module procedure :: dual_cos
  end interface

  intrinsic tan
  interface tan
    module procedure :: dual_tan
  end interface

contains

  subroutine dual_init(this, n, x, i)
    !! initialize dual number

    class(dual),       intent(out) :: this
    integer,           intent(in)  :: n
      !! number of base variables
    real,              intent(in)  :: x
      !! initial value
    integer, optional, intent(in)  :: i
      !! base variable index where 0 means no index (default: 0)

    ! set members
    this%n = n
    this%x = x

    ! init derivatives
    allocate (this%dx(n), source = 0.0)
    if (present(i)) then
      if (i .gt. 0) then
        ! set derivative wrt base variable i to 1
        this%dx(i) = 1
      end if
    end if
  end subroutine

  function dual_add_dual(x, y) result(r)
    !! add two dual numbers

    type(dual), intent(in) :: x
    type(dual), intent(in) :: y
    type(dual)             :: r

    r%n  = x%n
    r%x  = x%x + y%x
    r%dx = x%dx + y%dx
  end function

  function dual_add_real(x, y) result(r)
    !! add dual and real number

    type(dual), intent(in) :: x
    real,       intent(in) :: y
    type(dual)             :: r

    r%n  = x%n
    r%x  = x%x + y
    r%dx = x%dx
  end function

  function real_add_dual(x, y) result(r)
    !! add real and dual number

    real,       intent(in) :: x
    type(dual), intent(in) :: y
    type(dual)             :: r

    r = y + x
  end function

  function dual_neg(x) result(r)
    !! negate dual number

    type(dual), intent(in) :: x
    type(dual)             :: r

    r%n  = x%n
    r%x  = - x%x
    r%dx = - x%dx
  end function

  function dual_sub_dual(x, y) result(r)
    !! subtract two dual numbers

    type(dual), intent(in) :: x
    type(dual), intent(in) :: y
    type(dual)             :: r

    r = x + (-y)
  end function

  function dual_sub_real(x, y) result(r)
    !! subtract real from dual number

    type(dual), intent(in) :: x
    real,       intent(in) :: y
    type(dual)             :: r

    r = x + (-y)
  end function

  function real_sub_dual(x, y) result(r)
    !! subtract dual from real number

    real,       intent(in) :: x
    type(dual), intent(in) :: y
    type(dual)             :: r

    r = (-y) + x
  end function

  function dual_mul_dual(x, y) result(r)
    !! multiply two dual numbers

    type(dual), intent(in) :: x
    type(dual), intent(in) :: y
    type(dual)             :: r

    r%n  = x%n
    r%x  = x%x * y%x
    r%dx = x%dx * y%x + x%x * y%dx
  end function

  function dual_mul_real(x, y) result(r)
    !! multiply dual with real number

    type(dual), intent(in) :: x
    real,       intent(in) :: y
    type(dual)             :: r

    r%n  = x%n
    r%x  = x%x * y
    r%dx = x%dx * y
  end function

  function real_mul_dual(x, y) result(r)
    !! multiply real with dual number

    real,       intent(in) :: x
    type(dual), intent(in) :: y
    type(dual)             :: r

    r = y * x
  end function

  function dual_div_dual(x, y) result(r)
    !! divide two dual numbers

    type(dual), intent(in) :: x
    type(dual), intent(in) :: y
    type(dual)             :: r

    r%n = x%n
    r%x = x%x / y%x
    r%dx = x%dx / y%x - (x%x / y%x**2) * y%dx
  end function

  function dual_div_real(x, y) result(r)
    !! divide dual by real number

    type(dual), intent(in) :: x
    real,       intent(in) :: y
    type(dual)             :: r

    r%n  = x%n
    r%x  = x%x / y
    r%dx = x%dx / y
  end function

  function real_div_dual(x, y) result(r)
    !! divide real by dual number

    real,       intent(in) :: x
    type(dual), intent(in) :: y
    type(dual)             :: r

    r%n  = y%n
    r%x  = x / y%x
    r%dx = (- x / y%x**2) * y%dx
  end function

  function dual_pow_dual(x, y) result(r)
    !! raise a dual to the power of a dual number

    type(dual), intent(in) :: x
    type(dual), intent(in) :: y
    type(dual)             :: r

    r = exp(log(x) * y)
  end function

  function dual_pow_real(x, y) result(r)
    !! raise a dual to the power of a real number

    type(dual), intent(in) :: x
    real,       intent(in) :: y
    type(dual)             :: r

    r = exp(log(x) * y)
  end function

  function dual_pow_int(x, y) result(r)
    !! raise a dual number to the power of an integer

    type(dual), intent(in) :: x
    integer,    intent(in) :: y
    type(dual)             :: r

    ! local variables
    integer    :: p
    type(dual) :: r2

    ! for negative integers: compute inverse first
    if (y < 0) then
      r = 1.0 / x
      p = -y
    elseif (y > 0) then
      r = x
      p = y
    else
      ! x**0 = 1
      call r%init(x%n, 1.0)
      return
    end if

    ! compute low integer powers by repeated multiplication; high powers by real power function
    select case (p)
    case (2)
      r = r * r
    case (3)
      r2 = r * r
      r  = r2 * r
    case (4)
      r2 = r * r
      r  = r2 * r2
    case default
      r = r ** real(p)
    end select
  end function

  function real_pow_dual(x, y) result(r)
    !! raise real to the power of a dual number

    real,       intent(in) :: x
    type(dual), intent(in) :: y
    type(dual)             :: r

    r = exp(log(x) * y)
  end function

  function dual_abs(x) result(r)
    !! get absolute value of dual number

    type(dual), intent(in) :: x
    type(dual)             :: r

    if (x%x .ge. 0) then
      r = x
    else
      r = - x
    end if
  end function

  function dual_sqrt(x) result(r)
    !! compute square root of dual number

    type(dual), intent(in) :: x
    type(dual)             :: r

    r%n  = x%n
    r%x  = sqrt(x%x)
    r%dx = x%dx / (2 * r%x)
  end function

  function dual_exp(x) result(r)
    !! compute exponential function of dual number

    type(dual), intent(in) :: x
    type(dual)             :: r

    r%n  = x%n
    r%x  = exp(x%x)
    r%dx = r%x * x%dx
  end function

  function dual_log(x) result(r)
    !! compute natural logarithm of dual number

    type(dual), intent(in) :: x
    type(dual)             :: r

    r%n  = x%n
    r%x  = log(x%x)
    r%dx = x%dx / x%x
  end function

  function dual_sin(x) result(r)
    !! compute sine of dual number

    type(dual), intent(in) :: x
    type(dual)             :: r

    r%n  = x%n
    r%x  = sin(x%x)
    r%dx = cos(x%x) * x%dx
  end function

  function dual_cos(x) result(r)
    !! compute cosine of dual number

    type(dual), intent(in) :: x
    type(dual)             :: r

    r%n  = x%n
    r%x  =   cos(x%x)
    r%dx = - sin(x%x) * x%dx
  end function

  function dual_tan(x) result(r)
    !! compute tangent of dual number

    type(dual), intent(in) :: x
    type(dual)             :: r

    r%n  = x%n
    r%x  = tan(x%x)
    r%dx = x%dx / (cos(x%x)**2)
  end function

end module