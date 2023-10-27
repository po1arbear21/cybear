subroutine dual_{}N{}_init(this, x, i)
  !! initialize dual number (optionally as a base variable)
  class(dual_{}N),   intent(out) :: this
  real,              intent(in)  :: x
    !! value
  integer, optional, intent(in)  :: i
    !! optional: provide base variable index (default: 0)

  this%x  = x
  this%dx = 0
  if (present(i)) this%dx(i) = 1
end subroutine

elemental function dual_add_dual_{}N{}(x, y) result(r)
  !! add two dual numbers
  type(dual_{}N), intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x  + y%x
  r%dx = x%dx + y%dx
end function

elemental function dual_add_real_{}N{}(x, y) result(r)
  !! add dual and real number
  type(dual_{}N), intent(in) :: x
  real,           intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x  + y
  r%dx = x%dx
end function

elemental function dual_add_int_{}N{}(x, y) result(r)
  !! add dual and integer number
  type(dual_{}N), intent(in) :: x
  integer,        intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x  + y
  r%dx = x%dx
end function

elemental function real_add_dual_{}N{}(x, y) result(r)
  !! add real and dual number
  real,           intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x + y%x
  r%dx =     y%dx
end function

elemental function int_add_dual_{}N{}(x, y) result(r)
  !! add integer and dual number
  integer,        intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x + y%x
  r%dx =     y%dx
end function

elemental function dual_neg_{}N{}(x) result(r)
  !! negate dual number
  type(dual_{}N), intent(in) :: x
  type(dual_{}N)             :: r

  r%x  = - x%x
  r%dx = - x%dx
end function

elemental function dual_sub_dual_{}N{}(x, y) result(r)
  !! subtract two dual numbers
  type(dual_{}N), intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x  - y%x
  r%dx = x%dx - y%dx
end function

elemental function dual_sub_real_{}N{}(x, y) result(r)
  !! subtract real from dual number
  type(dual_{}N), intent(in) :: x
  real,           intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x  - y
  r%dx = x%dx
end function

elemental function dual_sub_int_{}N{}(x, y) result(r)
  !! subtract integer from dual number
  type(dual_{}N), intent(in) :: x
  integer,        intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x  - y
  r%dx = x%dx
end function

elemental function real_sub_dual_{}N{}(x, y) result(r)
  !! subtract dual from real number
  real,           intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x - y%x
  r%dx =   - y%dx
end function

elemental function int_sub_dual_{}N{}(x, y) result(r)
  !! subtract integer from real number
  integer,        intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x - y%x
  r%dx =   - y%dx
end function

elemental function dual_mul_dual_{}N{}(x, y) result(r)
  !! multiply two dual numbers
  type(dual_{}N), intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x * y%x
  r%dx = x%dx * y%x + x%x * y%dx
end function

elemental function dual_mul_real_{}N{}(x, y) result(r)
  !! multiply dual with real number
  type(dual_{}N), intent(in) :: x
  real,           intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x  * y
  r%dx = x%dx * y
end function

elemental function dual_mul_int_{}N{}(x, y) result(r)
  !! multiply dual with integer number
  type(dual_{}N), intent(in) :: x
  integer,        intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x  * y
  r%dx = x%dx * y
end function

elemental function real_mul_dual_{}N{}(x, y) result(r)
  !! multiply real with dual number
  real,           intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x * y%x
  r%dx = x * y%dx
end function

elemental function int_mul_dual_{}N{}(x, y) result(r)
  !! multiply integer with dual number
  integer,        intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x * y%x
  r%dx = x * y%dx
end function

elemental function dual_div_dual_{}N{}(x, y) result(r)
  !! divide two dual numbers
  type(dual_{}N), intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x / y%x
  r%dx = (x%dx - (x%x / y%x) * y%dx) / y%x
end function

elemental function dual_div_real_{}N{}(x, y) result(r)
  !! divide dual by real number
  type(dual_{}N), intent(in) :: x
  real,           intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x  / y
  r%dx = x%dx / y
end function

elemental function dual_div_int_{}N{}(x, y) result(r)
  !! divide dual by integer number
  type(dual_{}N), intent(in) :: x
  integer,        intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x  / y
  r%dx = x%dx / y
end function

elemental function real_div_dual_{}N{}(x, y) result(r)
  !! divide real by dual number
  real,           intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x / y%x
  r%dx = - (r%x * y%dx) / y%x
end function

elemental function int_div_dual_{}N{}(x, y) result(r)
  !! divide integer by dual number
  integer,        intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x / y%x
  r%dx = - (r%x * y%dx) / y%x
end function

elemental function dual_pow_dual_{}N{}(x, y) result(r)
  !! raise a dual to the power of a dual number
  type(dual_{}N), intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x**y%x
  r%dx = y%x * x%x ** (y%x-1) * x%dx

  if ((x%x /= 0) .or. (y%x <= 0)) then
    r%dx = r%dx + r%x * log(x%x) * y%dx
  end if
end function

elemental function dual_pow_real_{}N{}(x, y) result(r)
  !! raise a dual to the power of a real number
  type(dual_{}N), intent(in) :: x
  real,           intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x%x**y
  r%dx = y * x%x**(y-1) * x%dx
end function

elemental function dual_pow_int_{}N{}(x, y) result(r)
  !! raise a dual number to the power of an integer
  type(dual_{}N), intent(in) :: x
  integer,        intent(in) :: y
  type(dual_{}N)             :: r

  integer         :: p
  type(dual_{}N) :: r2

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

elemental function real_pow_dual_{}N{}(x, y) result(r)
  !! raise real to the power of a dual number
  real,           intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r%x  = x**y%x
  if ((x == 0) .and. (y%x > 0))then
    r%dx = 0 ! avoid log(0)
  else
    r%dx = r%x * log(x) * y%dx
  end if
end function

elemental function int_pow_dual_{}N{}(x, y) result(r)
  !! raise integer to the power of a dual number
  integer,        intent(in) :: x
  type(dual_{}N), intent(in) :: y
  type(dual_{}N)             :: r

  r = real(x) ** y
end function

function dual_dot_dual_{}N{}(x, y) result(r)
  !! dot product of two dual vectors
  type(dual_{}N), intent(in) :: x(:)
  type(dual_{}N), intent(in) :: y(:)
  type(dual_{}N)             :: r

  integer :: i

  m4_assert(size(x) == size(y))

  call r%init(0.0)
  do i = 1, size(x)
    r = r + x(i) * y(i)
  end do
end function

function real_dot_dual_{}N{}(x, y) result(r)
  !! dot product of real+dual vectors
  real,           intent(in) :: x(:)
  type(dual_{}N), intent(in) :: y(:)
  type(dual_{}N)             :: r

  integer :: i

  m4_assert(size(x) == size(y))

  call r%init(0.0)
  do i = 1, size(x)
    r = r + x(i) * y(i)
  end do
end function

function dual_dot_real_{}N{}(x, y) result(r)
  !! dot product of dual+real vectors
  type(dual_{}N), intent(in) :: x(:)
  real,           intent(in) :: y(:)
  type(dual_{}N)             :: r

  r = dot_product(y, x)
end function

elemental function dual_abs_{}N{}(x) result(r)
  !! get absolute value of dual number
  type(dual_{}N), intent(in) :: x
  type(dual_{}N)             :: r

  if (x%x >= 0) then
    r =  x
  else
    r = -x
  end if
end function

function dual_sum_{}N{}(x) result(r)
  !! get absolute value of dual number
  type(dual_{}N), intent(in) :: x(:)
  type(dual_{}N)             :: r

  integer :: i

  call r%init(0.0)
  do i = 1, size(x)
    r = r + x(i)
  end do
end function

elemental function dual_sqrt_{}N{}(x) result(r)
  !! compute square root of dual number
  type(dual_{}N), intent(in) :: x
  type(dual_{}N)             :: r

  r%x  = sqrt(x%x)
  r%dx = x%dx / (2 * r%x)
end function

elemental function dual_exp_{}N{}(x) result(r)
  !! compute exponential function of dual number
  type(dual_{}N), intent(in) :: x
  type(dual_{}N)             :: r

  r%x  = exp(x%x)
  r%dx = r%x * x%dx
end function

elemental function dual_expm1_{}N{}(x) result(r)
  !! compute exp(x) - 1
  type(dual_{}N), intent(in) :: x
  type(dual_{}N)             :: r

  r%x  = expm1(x%x)
  r%dx = exp(x%x) * x%dx
end function

elemental function dual_log_{}N{}(x) result(r)
  !! compute natural logarithm of dual number
  type(dual_{}N), intent(in) :: x
  type(dual_{}N)             :: r

  r%x  = log(x%x)
  r%dx = x%dx / x%x
end function

elemental function dual_log1p_{}N{}(x) result(r)
  !! compute log(x + 1)
  type(dual_{}N), intent(in) :: x
  type(dual_{}N)             :: r

  r%x  = log1p(x%x)
  r%dx = x%dx / (1.0 + x%x)
end function

elemental function dual_sin_{}N{}(x) result(r)
  !! compute sine of dual number
  type(dual_{}N), intent(in) :: x
  type(dual_{}N)             :: r

  r%x  = sin(x%x)
  r%dx = cos(x%x) * x%dx
end function

elemental function dual_cos_{}N{}(x) result(r)
  !! compute cosine of dual number
  type(dual_{}N), intent(in) :: x
  type(dual_{}N)             :: r

  r%x  = cos(x%x)
  r%dx = - sin(x%x) * x%dx
end function

elemental function dual_tan_{}N{}(x) result(r)
  !! compute tangent of dual number
  type(dual_{}N), intent(in) :: x
  type(dual_{}N)             :: r

  r%x  = tan(x%x)
  r%dx = x%dx / (cos(x%x)**2)
end function

m4_undefine({N})
