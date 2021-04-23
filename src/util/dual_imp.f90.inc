#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)

#define DUAL_TYPE            CONCAT(dual,N)
#define DUAL_INIT            CONCAT(dual_init,N)
#define DUAL_ADD_DUAL        CONCAT(dual_add_dual,N)
#define DUAL_ADD_REAL        CONCAT(dual_add_real,N)
#define DUAL_ADD_INT         CONCAT(dual_add_int,N)
#define REAL_ADD_DUAL        CONCAT(real_add_dual,N)
#define INT_ADD_DUAL         CONCAT(int_add_dual,N)
#define DUAL_NEG             CONCAT(dual_neg,N)
#define DUAL_SUB_DUAL        CONCAT(dual_sub_dual,N)
#define DUAL_SUB_REAL        CONCAT(dual_sub_real,N)
#define DUAL_SUB_INT         CONCAT(dual_sub_int,N)
#define REAL_SUB_DUAL        CONCAT(real_sub_dual,N)
#define INT_SUB_DUAL         CONCAT(int_sub_dual,N)
#define DUAL_MUL_DUAL        CONCAT(dual_mul_dual,N)
#define DUAL_MUL_REAL        CONCAT(dual_mul_real,N)
#define DUAL_MUL_INT         CONCAT(dual_mul_int,N)
#define REAL_MUL_DUAL        CONCAT(real_mul_dual,N)
#define INT_MUL_DUAL         CONCAT(int_mul_dual,N)
#define DUAL_DIV_DUAL        CONCAT(dual_div_dual,N)
#define DUAL_DIV_REAL        CONCAT(dual_div_real,N)
#define DUAL_DIV_INT         CONCAT(dual_div_int,N)
#define REAL_DIV_DUAL        CONCAT(real_div_dual,N)
#define INT_DIV_DUAL         CONCAT(int_div_dual,N)
#define DUAL_POW_DUAL        CONCAT(dual_pow_dual,N)
#define DUAL_POW_REAL        CONCAT(dual_pow_real,N)
#define DUAL_POW_INT         CONCAT(dual_pow_int,N)
#define REAL_POW_DUAL        CONCAT(real_pow_dual,N)
#define INT_POW_DUAL         CONCAT(int_pow_dual,N)
#define DUAL_DOT_DUAL        CONCAT(dual_dot_dual,N)
#define REAL_DOT_DUAL        CONCAT(real_dot_dual,N)
#define DUAL_DOT_REAL        CONCAT(dual_dot_real,N)
#define DUAL_ABS             CONCAT(dual_abs,N)
#define DUAL_SQRT            CONCAT(dual_sqrt,N)
#define DUAL_EXP             CONCAT(dual_exp,N)
#define DUAL_LOG             CONCAT(dual_log,N)
#define DUAL_SIN             CONCAT(dual_sin,N)
#define DUAL_COS             CONCAT(dual_cos,N)
#define DUAL_TAN             CONCAT(dual_tan,N)

subroutine DUAL_INIT(this, x, i)
  !! initialize dual number (optionally as a base variable)
  class(DUAL_TYPE),  intent(out) :: this
  real,              intent(in)  :: x
    !! value
  integer, optional, intent(in)  :: i
    !! optional: provide base variable index (default: 0)

  this%x  = x
  this%dx = 0
  if (present(i)) this%dx(i) = 1
end subroutine

elemental function DUAL_ADD_DUAL(x, y) result(r)
  !! add two dual numbers
  type(DUAL_TYPE), intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x  + y%x
  r%dx = x%dx + y%dx
end function

elemental function DUAL_ADD_REAL(x, y) result(r)
  !! add dual and real number
  type(DUAL_TYPE), intent(in) :: x
  real,            intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x  + y
  r%dx = x%dx
end function

elemental function DUAL_ADD_INT(x, y) result(r)
  !! add dual and integer number
  type(DUAL_TYPE), intent(in) :: x
  integer,         intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x  + y
  r%dx = x%dx
end function

elemental function REAL_ADD_DUAL(x, y) result(r)
  !! add real and dual number
  real,            intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x + y%x
  r%dx =     y%dx
end function

elemental function INT_ADD_DUAL(x, y) result(r)
  !! add integer and dual number
  integer,         intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x + y%x
  r%dx =     y%dx
end function

elemental function DUAL_NEG(x) result(r)
  !! negate dual number
  type(DUAL_TYPE), intent(in) :: x
  type(DUAL_TYPE)             :: r

  r%x  = - x%x
  r%dx = - x%dx
end function

elemental function DUAL_SUB_DUAL(x, y) result(r)
  !! subtract two dual numbers
  type(DUAL_TYPE), intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x  - y%x
  r%dx = x%dx - y%dx
end function

elemental function DUAL_SUB_REAL(x, y) result(r)
  !! subtract real from dual number
  type(DUAL_TYPE), intent(in) :: x
  real,            intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x  - y
  r%dx = x%dx
end function

elemental function DUAL_SUB_INT(x, y) result(r)
  !! subtract integer from dual number
  type(DUAL_TYPE), intent(in) :: x
  integer,         intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x  - y
  r%dx = x%dx
end function

elemental function REAL_SUB_DUAL(x, y) result(r)
  !! subtract dual from real number
  real,            intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x - y%x
  r%dx =   - y%dx
end function

elemental function INT_SUB_DUAL(x, y) result(r)
  !! subtract integer from real number
  integer,         intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x - y%x
  r%dx =   - y%dx
end function

elemental function DUAL_MUL_DUAL(x, y) result(r)
  !! multiply two dual numbers
  type(DUAL_TYPE), intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x * y%x
  r%dx = x%dx * y%x + x%x * y%dx
end function

elemental function DUAL_MUL_REAL(x, y) result(r)
  !! multiply dual with real number
  type(DUAL_TYPE), intent(in) :: x
  real,            intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x  * y
  r%dx = x%dx * y
end function

elemental function DUAL_MUL_INT(x, y) result(r)
  !! multiply dual with integer number
  type(DUAL_TYPE), intent(in) :: x
  integer,         intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x  * y
  r%dx = x%dx * y
end function

elemental function REAL_MUL_DUAL(x, y) result(r)
  !! multiply real with dual number
  real,            intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x * y%x
  r%dx = x * y%dx
end function

elemental function INT_MUL_DUAL(x, y) result(r)
  !! multiply integer with dual number
  integer,         intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x * y%x
  r%dx = x * y%dx
end function

elemental function DUAL_DIV_DUAL(x, y) result(r)
  !! divide two dual numbers
  type(DUAL_TYPE), intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x / y%x
  r%dx = x%dx / y%x - (x%x / y%x**2) * y%dx
end function

elemental function DUAL_DIV_REAL(x, y) result(r)
  !! divide dual by real number
  type(DUAL_TYPE), intent(in) :: x
  real,            intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x  / y
  r%dx = x%dx / y
end function

elemental function DUAL_DIV_INT(x, y) result(r)
  !! divide dual by integer number
  type(DUAL_TYPE), intent(in) :: x
  integer,         intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x  / y
  r%dx = x%dx / y
end function

elemental function REAL_DIV_DUAL(x, y) result(r)
  !! divide real by dual number
  real,            intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x / y%x
  r%dx = (- x / y%x**2) * y%dx
end function

elemental function INT_DIV_DUAL(x, y) result(r)
  !! divide integer by dual number
  integer,         intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x / y%x
  r%dx = (- x / y%x**2) * y%dx
end function

elemental function DUAL_POW_DUAL(x, y) result(r)
  !! raise a dual to the power of a dual number
  type(DUAL_TYPE), intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x**y%x
  r%dx = y%x * x%x ** (y%x-1) * x%dx

  if ((x%x /= 0) .or. (y%x <= 0)) then
    r%dx = r%dx + r%x * log(x%x) * y%dx
  end if
end function

elemental function DUAL_POW_REAL(x, y) result(r)
  !! raise a dual to the power of a real number
  type(DUAL_TYPE), intent(in) :: x
  real,            intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x%x**y
  r%dx = y * x%x**(y-1) * x%dx
end function

elemental function DUAL_POW_INT(x, y) result(r)
  !! raise a dual number to the power of an integer
  type(DUAL_TYPE), intent(in) :: x
  integer,         intent(in) :: y
  type(DUAL_TYPE)             :: r

  integer         :: p
  type(DUAL_TYPE) :: r2

  ! for negative integers: compute inverse first
  if (y < 0) then
    r = 1.0 / x
    p = - y
  elseif (y > 0) then
    r = x
    p = y
  else
    r%x  = 1.0
    r%dx = 0.0
    return
  end if

  ! compute low integer powers by repeated multiplication; high powers by real power function
  ! FIXME: is this actually faster than just the default case?
  select case (p)
    case (2)
      r = r * r
    case (3)
      r = (r * r) * r
    case (4)
      r2 = r * r
      r  = r2 * r2
    case default
      r = r ** real(p)
  end select
end function

elemental function REAL_POW_DUAL(x, y) result(r)
  !! raise real to the power of a dual number
  real,            intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r%x  = x**y%x
  if ((x == 0) .and. (y%x > 0))then
    r%dx = 0 ! avoid log(0)
  else
    r%dx = r%x * log(x) * y%dx
  end if
end function

elemental function INT_POW_DUAL(x, y) result(r)
  !! raise integer to the power of a dual number
  integer,         intent(in) :: x
  type(DUAL_TYPE), intent(in) :: y
  type(DUAL_TYPE)             :: r

  r = real(x) ** y
end function

function DUAL_DOT_DUAL(x, y) result(r)
  !! dot product of two dual vectors
  type(DUAL_TYPE), intent(in) :: x(:)
  type(DUAL_TYPE), intent(in) :: y(:)
  type(DUAL_TYPE)             :: r

  integer :: i

  ASSERT(size(x) == size(y))

  r%x  = 0
  r%dx = 0
  do i = 1, size(x)
    r = r + x(i) * y(i)
  end do
end function

function REAL_DOT_DUAL(x, y) result(r)
  !! dot product of real+dual vectors
  real,            intent(in) :: x(:)
  type(DUAL_TYPE), intent(in) :: y(:)
  type(DUAL_TYPE)             :: r

  integer         :: i
  type(DUAL_TYPE) :: xx(size(x))

  do i = 1, size(x)
    call xx(i)%init(x(i))
  end do
  r = dot_product(xx, y)
end function

function DUAL_DOT_REAL(x, y) result(r)
  !! dot product of dual+real vectors
  type(DUAL_TYPE), intent(in) :: x(:)
  real,            intent(in) :: y(:)
  type(DUAL_TYPE)             :: r

  integer         :: i
  type(DUAL_TYPE) :: yy(size(y))

  do i = 1, size(y)
    call yy(i)%init(y(i))
  end do
  r = dot_product(x, yy)
end function

elemental function DUAL_ABS(x) result(r)
  !! get absolute value of dual number
  type(DUAL_TYPE), intent(in) :: x
  type(DUAL_TYPE)             :: r

  if (x%x >= 0) then
    r =  x
  else
    r = -x
  end if
end function

elemental function DUAL_SQRT(x) result(r)
  !! compute square root of dual number
  type(DUAL_TYPE), intent(in) :: x
  type(DUAL_TYPE)             :: r

  r%x  = sqrt(x%x)
  r%dx = x%dx / (2 * r%x)
end function

elemental function DUAL_EXP(x) result(r)
  !! compute exponential function of dual number
  type(DUAL_TYPE), intent(in) :: x
  type(DUAL_TYPE)             :: r

  r%x  = exp(x%x)
  r%dx = r%x * x%dx
end function

elemental function DUAL_LOG(x) result(r)
  !! compute natural logarithm of dual number
  type(DUAL_TYPE), intent(in) :: x
  type(DUAL_TYPE)             :: r

  r%x  = log(x%x)
  r%dx = x%dx / x%x
end function

elemental function DUAL_SIN(x) result(r)
  !! compute sine of dual number
  type(DUAL_TYPE), intent(in) :: x
  type(DUAL_TYPE)             :: r

  r%x  = sin(x%x)
  r%dx = cos(x%x) * x%dx
end function

elemental function DUAL_COS(x) result(r)
  !! compute cosine of dual number
  type(DUAL_TYPE), intent(in) :: x
  type(DUAL_TYPE)             :: r

  r%x  = cos(x%x)
  r%dx = - sin(x%x) * x%dx
end function

elemental function DUAL_TAN(x) result(r)
  !! compute tangent of dual number
  type(DUAL_TYPE), intent(in) :: x
  type(DUAL_TYPE)             :: r

  r%x  = tan(x%x)
  r%dx = x%dx / (cos(x%x)**2)
end function

#undef N

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef DUAL_TYPE
#undef DUAL_INIT
#undef DUAL_ADD_DUAL
#undef DUAL_ADD_REAL
#undef DUAL_ADD_INT
#undef REAL_ADD_DUAL
#undef INT_ADD_DUAL
#undef DUAL_NEG
#undef DUAL_SUB_DUAL
#undef DUAL_SUB_REAL
#undef DUAL_SUB_INT
#undef REAL_SUB_DUAL
#undef INT_SUB_DUAL
#undef DUAL_MUL_DUAL
#undef DUAL_MUL_REAL
#undef DUAL_MUL_INT
#undef REAL_MUL_DUAL
#undef INT_MUL_DUAL
#undef DUAL_DIV_DUAL
#undef DUAL_DIV_REAL
#undef DUAL_DIV_INT
#undef REAL_DIV_DUAL
#undef INT_DIV_DUAL
#undef DUAL_POW_DUAL
#undef DUAL_POW_REAL
#undef DUAL_POW_INT
#undef REAL_POW_DUAL
#undef INT_POW_DUAL
#undef DUAL_DOT_DUAL
#undef REAL_DOT_DUAL
#undef DUAL_DOT_REAL
#undef DUAL_ABS
#undef DUAL_SQRT
#undef DUAL_EXP
#undef DUAL_LOG
#undef DUAL_SIN
#undef DUAL_COS
#undef DUAL_TAN
