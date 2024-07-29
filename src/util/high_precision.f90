m4_include(macro.f90.inc)

module high_precision_m
  !! algorithms taken from:
  !!   "Accurate sum and dot prodcut" by Ogita, Rump and Oishi
  !!   "High precision evaluation of nonlinear functions" by Rump
  use error_m,         only: assert_failed, program_error
  use ieee_arithmetic

  implicit none

  private
  public hp_real
  public hp_to_real, real_to_hp
  public operator(+), operator(-), operator(*), operator(/), operator(**)
  public abs, sqrt, exp, expm1, log, log1p, sin, cos, tan, sinh, cosh, tanh, atanh, ber
  public hp_sum, hp_dot
  public TwoSum, TwoProduct, TwoDivision

  type hp_real
    !! represents high precision value x + y

    real :: x = 0.0
      !! principal value
    real :: y = 0.0
      !! correction value
  end type

  interface operator (+)
    module procedure :: hp_real_add_hh
    module procedure :: hp_real_add_hr
    module procedure :: hp_real_add_rh
    module procedure :: hp_real_add_hi
    module procedure :: hp_real_add_ih
  end interface

  interface operator (-)
    module procedure :: hp_real_neg
    module procedure :: hp_real_sub_hh
    module procedure :: hp_real_sub_hr
    module procedure :: hp_real_sub_rh
    module procedure :: hp_real_sub_hi
    module procedure :: hp_real_sub_ih
  end interface

  interface operator (*)
    module procedure :: hp_real_mul_hh
    module procedure :: hp_real_mul_hr
    module procedure :: hp_real_mul_rh
    module procedure :: hp_real_mul_hi
    module procedure :: hp_real_mul_ih
  end interface

  interface operator (/)
    module procedure :: hp_real_div_hh
    module procedure :: hp_real_div_hr
    module procedure :: hp_real_div_rh
    module procedure :: hp_real_div_hi
    module procedure :: hp_real_div_ih
  end interface

  interface operator (**)
    module procedure :: hp_real_pow_hh
    module procedure :: hp_real_pow_hr
    module procedure :: hp_real_pow_rh
    module procedure :: hp_real_pow_hi
    module procedure :: hp_real_pow_ih
  end interface

  interface abs
    module procedure :: hp_abs
  end interface

  interface sqrt
    module procedure :: hp_sqrt
  end interface

  interface exp
    module procedure :: hp_exp
  end interface

  interface expm1
    module procedure :: hp_expm1
  end interface

  interface log
    module procedure :: hp_log
  end interface

  interface log1p
    module procedure :: hp_log1p
  end interface

  interface sin
    module procedure :: hp_sin
  end interface

  interface cos
    module procedure :: hp_cos
  end interface

  interface tan
    module procedure :: hp_tan
  end interface

  interface sinh
    module procedure :: hp_sinh
  end interface

  interface cosh
    module procedure :: hp_cosh
  end interface

  interface tanh
    module procedure :: hp_tanh
  end interface

  interface atanh
    module procedure :: hp_atanh
  end interface

  interface ber
    module procedure :: hp_ber
  end interface

  interface hp_sum
    module procedure :: SumKvert
  end interface

  interface hp_dot
    module procedure :: DotK
  end interface

contains

  elemental function TwoSum(a, b) result(h)
    !! error free transformation of the sum of two floating point numbers
    real, intent(in) :: a
    real, intent(in) :: b
    type(hp_real)    :: h

    real             :: z

    ! principal
    h%x = a + b

    ! correction
    z   = h%x - a
    h%y = (a - (h%x - z)) + (b - z)
  end function

  elemental function Split(a) result(h)
    !! error free splitting of float into two parts
    real, intent(in) :: a
    type(hp_real)    :: h

    real            :: c
    real, parameter :: factor = 134217729.0      ! 2^27+1

    c   = factor * a
    h%x = c - (c - a)
    h%y = a - h%x
  end function

  elemental function SplitQuad(a) result(h)
    !! error free splitting of 128 bit float into two 64 bit floats
    real(kind=16), intent(in) :: a
    type(hp_real)             :: h

    real(kind=16)            :: c
    real(kind=16), parameter :: factor = 18014398509481985.0 ! 2^54+1

    c   = factor * a
    h%x = real(c - (c - a))
    h%y = real(a - h%x)
  end function

  elemental function TwoProduct(a, b) result(h)
    !! multiply two real numbers with high precision (error free product)
    real, intent(in) :: a
    real, intent(in) :: b
    type(hp_real)    :: h

    type(hp_real) :: h1, h2

    ! principal
    h%x = a * b

    ! correction
    h1 = Split(a)
    h2 = Split(b)
    associate(y => h%y, x=> h%x, a1 => h1%x, a2 => h1%y, b1 => h2%x, b2 => h2%y)
      y = a2*b2 - (((x-a1*b1) - a2*b1) - a1*b2)
    end associate
  end function

  elemental function TwoDivision(a, b) result(h)
    !! divide two real numbers with high precision (not error free)
    real, intent(in) :: a
    real, intent(in) :: b
    type(hp_real)    :: h

    real :: x0

    ! first approximation
    x0 = a / b

    ! residual
    h = a - TwoProduct(b, x0)

    ! Newton step
    h%x = h%x / b
    h%y = h%y / b
    h = x0 + h
  end function

  elemental function real_to_hp(r) result(h)
    !! convert real to high precision real
    real, intent(in) :: r
    type(hp_real)    :: h

    ! principal is r; correction zero
    h%x = r
    h%y = 0
  end function

  elemental function hp_to_real(h) result(r)
    !! convert high precision real to real
    type(hp_real), intent(in) :: h
    real                      :: r

    ! add principal and correction
    r = h%x + h%y
  end function

  elemental function hp_real_add_hh(h1, h2) result(h3)
    !! add two high precision reals
    type(hp_real), intent(in) :: h1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    real :: c

    h3   = TwoSum(h1%x, h2%x)
    c    = h3%y
    h3   = TwoSum(h3%x, h1%y)
    c    = c + h3%y
    h3   = TwoSum(h3%x, h2%y)
    h3%y = h3%y + c
  end function

  elemental function hp_real_add_hr(h1, r2) result(h3)
    !! add real to high precision real
    type(hp_real), intent(in) :: h1
    real,          intent(in) :: r2
    type(hp_real)             :: h3

    real :: c

    h3   = TwoSum(h1%x, r2)
    c    = h3%y
    h3   = TwoSum(h3%x, h1%y)
    h3%y = h3%y + c
  end function

  elemental function hp_real_add_rh(r1, h2) result(h3)
    !! add high precision real to real
    real,          intent(in) :: r1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = h2 + r1 ! calls hp_real_add_hr
  end function

  elemental function hp_real_add_hi(h1, i2) result(h3)
    !! add integer to high precision real
    type(hp_real), intent(in) :: h1
    integer,       intent(in) :: i2
    type(hp_real)             :: h3

    h3 = h1 + real(i2)
  end function

  elemental function hp_real_add_ih(i1, h2) result(h3)
    !! add high precision real to integer
    integer,       intent(in) :: i1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = real(i1) + h2
  end function

  elemental function hp_real_neg(h1) result(h2)
    !! negate high precision real
    type(hp_real), intent(in) :: h1
    type(hp_real)             :: h2

    ! negate principal and correction
    h2%x = - h1%x
    h2%y = - h1%y
  end function

  elemental function hp_real_sub_hh(h1, h2) result(h3)
    !! subtract two high precision reals
    type(hp_real), intent(in) :: h1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = h1 + (-h2) ! calls hp_real_neg -> hp_real_add_hh
  end function

  elemental function hp_real_sub_hr(h1, r2) result(h3)
    !! subtract real from high precision real
    type(hp_real), intent(in) :: h1
    real,          intent(in) :: r2
    type(hp_real)             :: h3

    h3 = h1 + (-r2)           ! calls hp_real_add_hr
  end function

  elemental function hp_real_sub_rh(r1, h2) result(h3)
    !! subtract high precision real from real
    real,          intent(in) :: r1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = r1 + (-h2) ! calls hp_real_neg, hp_real_add_rh
  end function

  elemental function hp_real_sub_hi(h1, i2) result(h3)
    !! subtract integer from high precision real
    type(hp_real), intent(in) :: h1
    integer,       intent(in) :: i2
    type(hp_real)             :: h3

    h3 = h1 - real(i2)
  end function

  elemental function hp_real_sub_ih(i1, h2) result(h3)
    !! subtract high precision real from integer
    integer,       intent(in) :: i1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = real(i1) - h2
  end function

  elemental function hp_real_mul_hh(h1, h2) result(h3)
    !! multiply two high precision reals
    type(hp_real), intent(in) :: h1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    real :: c

    h3   = TwoProduct(h1%x, h2%x)
    c    = h3%y
    h3   = TwoSum(h3%x, h1%x*h2%y)
    c    = c + h3%y
    h3   = TwoSum(h3%x, h1%y*h2%x)
    h3%y = h3%y + c
  end function

  elemental function hp_real_mul_hr(h1, r2) result(h3)
    !! multiply high precision real by real
    type(hp_real), intent(in) :: h1
    real,          intent(in) :: r2
    type(hp_real)             :: h3

    real :: c

    h3   = TwoProduct(h1%x, r2)
    c    = h3%y
    h3   = TwoSum(h3%x, h1%y*r2)
    h3%y = h3%y + c
  end function

  elemental function hp_real_mul_rh(r1, h2) result(h3)
    !! multiply real by high precision real
    real,          intent(in) :: r1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = h2 * r1
  end function

  elemental function hp_real_mul_hi(h1, i2) result(h3)
    !! multiply high precision real by integer
    type(hp_real), intent(in) :: h1
    integer,       intent(in) :: i2
    type(hp_real)             :: h3

    h3 = h1 * real(i2)
  end function

  elemental function hp_real_mul_ih(i1, h2) result(h3)
    !! multiply integer by high precision real
    integer,       intent(in) :: i1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = real(i1) * h2
  end function

  elemental function hp_real_div_hh(h1, h2) result(h3)
    !! divide two high precision reals
    type(hp_real), intent(in) :: h1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    real :: c

    ! first approximation
    c = h1%x / h2%x

    ! one step of newton iteration
    h3 = c - (h2 * c - h1) / h2%x
  end function

  elemental function hp_real_div_hr(h1, r2) result(h3)
    !! divide high precision real by real
    type(hp_real), intent(in) :: h1
    real,          intent(in) :: r2
    type(hp_real)             :: h3

    h3 = TwoDivision(h1%x, r2) + TwoDivision(h1%y, r2)
  end function

  elemental function hp_real_div_rh(r1, h2) result(h3)
    !! divide real by high precision real
    real,          intent(in) :: r1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    real :: c

    ! first approximation
    c = r1 / h2%x

    ! one step of newton iteration
    h3 = c - (h2 * c - r1) / h2%x
  end function

  elemental function hp_real_div_hi(h1, i2) result(h3)
    !! divide high precision real by integer
    type(hp_real), intent(in) :: h1
    integer,       intent(in) :: i2
    type(hp_real)             :: h3

    h3 = h1 / real(i2)
  end function

  elemental function hp_real_div_ih(i1, h2) result(h3)
    !! divide integer by high precision real
    integer,       intent(in) :: i1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = real(i1) / h2
  end function

  elemental function hp_real_pow_hh(h1, h2) result(h3)
    !! h3 = h1 ** h2
    type(hp_real), intent(in) :: h1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    if (hp_to_real(h2) == 0) then
      h3 = real_to_hp(1.0)
    elseif (hp_to_real(h1) == 0) then
      h3 = real_to_hp(0.0)
    else
      h3 = exp(log(h1) * h2)
    end if
  end function

  elemental function hp_real_pow_hr(h1, r2) result(h3)
    !! h3 = h1 ** r2
    type(hp_real), intent(in) :: h1
    real,          intent(in) :: r2
    type(hp_real)             :: h3

    if (r2 == 0) then
      h3 = real_to_hp(1.0)
    elseif (hp_to_real(h1) == 0) then
      h3 = real_to_hp(0.0)
    else
      h3 = exp(log(h1) * r2)
    end if
  end function

  elemental function hp_real_pow_rh(r1, h2) result(h3)
    !! h3 = r1 ** h2
    real,          intent(in) :: r1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    if (hp_to_real(h2) == 0) then
      h3 = real_to_hp(1.0)
    elseif (r1 == 0) then
      h3 = real_to_hp(0.0)
    else
      h3 = exp(log(real_to_hp(r1)) * h2)
    end if
  end function

  elemental function hp_real_pow_hi(h1, i2) result(h3)
    !! h3 = h1 ** i2
    type(hp_real), intent(in) :: h1
    integer,       intent(in) :: i2
    type(hp_real)             :: h3

    select case (i2)
    case (-3)
      h3 = 1.0 / (h1 * h1 * h1)
    case (-2)
      h3 = 1.0 / (h1 * h1)
    case (-1)
      h3 = 1.0 / h1
    case ( 0)
      h3 = real_to_hp(1.0)
    case (+1)
      h3 = h1
    case (+2)
      h3 = h1 * h1
    case (+3)
      h3 = h1 * h1 * h1
    case default
      h3 = (sign(1.0, hp_to_real(h1))**i2) * abs(h1) ** real(i2)
    end select
  end function

  elemental function hp_real_pow_ih(i1, h2) result(h3)
    !! h3 = i1 ** h2
    integer,       intent(in) :: i1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = real(i1) ** h2
  end function

  elemental function hp_abs(h) result(a)
    !! absolute value of high precision real
    type(hp_real), intent(in) :: h
    type(hp_real)             :: a

    if (hp_to_real(h) >= 0.0) then
      a = h
    else
      a = -h
    end if
  end function

  elemental function hp_sqrt(h) result(s)
    !! high precision square root
    type(hp_real), intent(in) :: h
    type(hp_real)             :: s

    real :: s0

    ! first approximation
    s0 = sqrt(hp_to_real(h))
    s = real_to_hp(s0)

    ! one step of newton iteration
    if (s0 > 0) then
      s = 0.5 * (s + h / s)
    end if
  end function

  elemental function hp_exp(h) result(e)
    !! high precision exponential function (simply uses quadruple precision)
    type(hp_real), intent(in) :: h
    type(hp_real)             :: e

    real(kind=16) :: tmp

    tmp = h%x
    tmp = tmp + h%y
    tmp = exp(tmp)

    e = SplitQuad(tmp)
  end function

  elemental function hp_expm1(h) result(e)
    !! high precision exp(h) - 1
    type(hp_real), intent(in) :: h
    type(hp_real)             :: e

    real(kind=16) :: tmp, etmp

    if (ieee_class(h%x) == IEEE_POSITIVE_INF) then
      e%x = h%x
      e%y = 0.0
      return
    end if
    if (ieee_class(h%y) == IEEE_POSITIVE_INF) then
      e%x = h%y
      e%y = 0.0
      return
    end if

    tmp = h%x
    tmp = tmp + h%y

    etmp = exp(tmp)

    if (etmp == 1.0) then
      e = h
    else if (etmp - 1.0 == -1.0) then
      e%x = -1.0
      e%y =  0.0
    else
      tmp = (etmp - 1.0) * tmp / log(etmp)
      e = SplitQuad(tmp)
    end if
  end function

  elemental function hp_log(h) result(l)
    !! high precision logarithm (simply uses quadruple precision)
    type(hp_real), intent(in) :: h
    type(hp_real)             :: l

    real(kind=16) :: tmp

    tmp = h%x
    tmp = tmp + h%y
    tmp = log(tmp)

    l = SplitQuad(tmp)
  end function

  elemental function hp_log1p(h) result(l)
    !! high precision log(1 + x)
    type(hp_real), intent(in) :: h
    type(hp_real)             :: l

    real(kind=16) :: tmp, u, d

    if (ieee_class(h%x) == IEEE_POSITIVE_INF) then
      l%x = h%x
      l%y = 0.0
      return
    end if
    if (ieee_class(h%y) == IEEE_POSITIVE_INF) then
      l%x = h%y
      l%y = 0.0
      return
    end if

    tmp = h%x
    tmp = tmp + h%y

    u = tmp + 1.0
    d = u - 1.0

    if (d == 0) then
      l = h
    else
      tmp = log(u) * tmp / d
      l = SplitQuad(tmp)
    end if
  end function

  elemental function hp_sin(h) result(s)
    !! high precision sine function (simply uses quadruple precision)
    type(hp_real), intent(in) :: h
    type(hp_real)             :: s

    real(kind=16) :: tmp

    tmp = h%x
    tmp = tmp + h%y
    tmp = sin(tmp)

    s = SplitQuad(tmp)
  end function

  elemental function hp_cos(h) result(c)
    !! high precision cosine function (simply uses quadruple precision)
    type(hp_real), intent(in) :: h
    type(hp_real)             :: c

    real(kind=16) :: tmp

    tmp = h%x
    tmp = tmp + h%y
    tmp = cos(tmp)

    c = SplitQuad(tmp)
  end function

  elemental function hp_tan(h) result(t)
    !! high precision tangent function (simply uses quadruple precision)
    type(hp_real), intent(in) :: h
    type(hp_real)             :: t

    real(kind=16) :: tmp

    tmp = h%x
    tmp = tmp + h%y
    tmp = tan(tmp)

    t = SplitQuad(tmp)
  end function

  elemental function hp_sinh(h) result(s)
    !! high precision hyperbolic sine function (simply uses quadruple precision)
    type(hp_real), intent(in) :: h
    type(hp_real)             :: s

    real(kind=16) :: tmp

    tmp = h%x
    tmp = tmp + h%y
    tmp = sinh(tmp)

    s = SplitQuad(tmp)
  end function

  elemental function hp_cosh(h) result(c)
    !! high precision hyperbolic cosine function (simply uses quadruple precision)
    type(hp_real), intent(in) :: h
    type(hp_real)             :: c

    real(kind=16) :: tmp

    tmp = h%x
    tmp = tmp + h%y
    tmp = cosh(tmp)

    c = SplitQuad(tmp)
  end function

  elemental function hp_tanh(h) result(t)
    !! high precision hyperbolic tangent function (simply uses quadruple precision)
    type(hp_real), intent(in) :: h
    type(hp_real)             :: t

    real(kind=16) :: tmp

    tmp = h%x
    tmp = tmp + h%y
    tmp = tanh(tmp)

    t = SplitQuad(tmp)
  end function

  elemental function hp_atanh(h) result(t)
    !! high precision inverse hyperbolic tangent function (simply uses quadruple precision)
    type(hp_real), intent(in) :: h
    type(hp_real)             :: t

    real(kind=16) :: tmp

    tmp = h%x
    tmp = tmp + h%y
    tmp = atanh(tmp)

    t = SplitQuad(tmp)
  end function

  elemental function hp_ber(h) result(b)
    !! high precision bernoulli function (simply uses quadruple precision)
    type(hp_real), intent(in) :: h
    type(hp_real)             :: b

    real(kind=16) :: tmp

    tmp = h%x
    tmp = tmp + h%y

    if (abs(tmp) > 1e-8) then
      tmp = 0.5 * tmp * exp(-0.5 * tmp) / sinh(0.5 * tmp)
    else
      tmp = 1.0 + tmp * (-0.5 + tmp / 12.0)
    end if

    b = SplitQuad(tmp)
  end function

  function SumKvert(p, K) result(res)
    !! high precision K-fold summation.
    real              , intent(in)  :: p(:)
    integer, optional , intent(in)  :: K
      !! K-fold precision.
      !! default: 2
    real                            :: res

    integer           :: n, i, j, k_, KK
    real              :: s, alp
    real, allocatable :: q(:)
    type(hp_real)     :: h_tmp

    m4_assert(size(p) >= 0)

    ! optional arg
    KK = 2
    if (present(K)) KK = K

    ! handle K special case
    if      (KK <  1) then
      call program_error('K should be greater than 0')
    else if (KK == 1) then
      res = sum(p)
      return
    end if

    ! handle p special case
    if (size(p) < 2) then
      res = sum(p)
      return
    end if

    n  = size(p)
    KK = min(KK, n)
    allocate (q(KK-1))

    do i = 1, KK-1
      s = p(i)
      do k_ = 1, i-1
        ! [qk, s] <- TowSum(qk, s)
        h_tmp = TwoSum(q(k_), s)
        q(k_) = h_tmp%x
        s     = h_tmp%y
      end do
      q(i) = s
    end do

    s = 0.0
    do i = KK, n
      alp = p(i)
      do k_ = 1, KK-1
        ! [qk, alpha] <- TowSum(qk, alpha)
        h_tmp = TwoSum(q(k_), alp)
        q(k_) = h_tmp%x
        alp   = h_tmp%y
      end do
      s = s + alp
    end do

    do j = 1, KK-2
      alp = q(j)
      do k_ = j+1, KK-1
        ! [qk, alpha] <- TowSum(qk, alpha)
        h_tmp = TwoSum(q(k_), alp)
        q(k_) = h_tmp%x
        alp   = h_tmp%y
      end do
      s = s + alp
    end do

    res = s + q(KK-1)
  end function

  function DotK(x, y, K) result(res)
    !! high precision K-fold dot product.
    real              , intent(in)  :: x(:)
    real              , intent(in)  :: y(:)
    integer, optional , intent(in)  :: K
      !! K-fold precision.
      !! default: 3
    real                            :: res

    integer           :: n, i, K_
    real              :: p, h
    real, allocatable :: r(:)
    type(hp_real)     :: hp_tmp

    K_ = 3
    if (present(K)) K_ = K

    if (K_ <  3) call program_error('K must be greater than 2')

    n = size(x)
    allocate(r(2*n))

    ! [p, r1] = TwoProduct(x1, y1)
    hp_tmp  = TwoProduct(x(1), y(1))
    p       = hp_tmp%x
    r(1)    = hp_tmp%y

    do i = 2, n
      ! [h, ri] = TwoProduct(xi, yi)
      hp_tmp  = TwoProduct(x(i), y(i))
      h       = hp_tmp%x
      r(i)    = hp_tmp%y

      ! [p, r_{n+i-1}] = TwoSum(p, h)
      hp_tmp    = TwoSum(p, h)
      p         = hp_tmp%x
      r(n+i-1)  = hp_tmp%y
    end do
    r(2*n) = p
    res = SumKvert(r, K_-1)
  end function

end module
