#include "macro.f90.inc"

module math_m
  use error_m
  use ieee_arithmetic
  implicit none

  private
  public :: PI
  public :: cross_product, cross_product_2d
  public :: heaviside
  public :: isinf
  public :: ber, dberdx
  public :: phi1, dphi1dx, phi2, dphi2dx
  public :: expm1, log1p
  public :: linspace, logspace
  public :: eye_int, eye_real
  public :: norm_inf
  public :: check_lin_dep

  real, parameter :: PI = 3.141592653589793238462643

  interface expm1
    module procedure :: m_expm1
  end interface
  interface log1p
    module procedure :: m_log1p
  end interface

contains

  pure function cross_product(a, b) result(axb)
    real, intent(in) :: a(3), b(3)
    real             :: axb(3)

    axb = [a(2)*b(3) - a(3)*b(2), &
      &    a(3)*b(1) - a(1)*b(3), &
      &    a(1)*b(2) - a(2)*b(1)  ]
  end function

  pure real function cross_product_2d(a, b)
    real, intent(in) :: a(2), b(2)

    cross_product_2d = a(1)*b(2) - a(2)*b(1)
  end function

  elemental function heaviside(x) result(h)
    !! heaviside step function
    real, intent(in) :: x
    real             :: h

    h = 0.5 * (sign(1.0, x) + 1.0)
  end function

  elemental function isinf(x) result (r)
    !! return true if argument is positive of negative infinity
    real, intent(in) :: x
    logical          :: r

    r = (ieee_class(x) == ieee_positive_inf) .or. (ieee_class(x) == ieee_negative_inf)
  end function

  elemental function ber(x) result(b)
    !! bernoulli function
    real, intent(in) :: x
    real             :: b

    if (abs(x) > 1e-6) then
      b = 0.5 * x * exp(-0.5 * x) / sinh(0.5 * x)
    else
      b = 1.0 + x * (-0.5 + x / 12.0)
    end if
  end function

  elemental function dberdx(x) result(dbdx)
    !! derivative of bernoulli function
    real, intent(in) :: x
    real             :: dbdx

    if (abs(x) > 1e-6) then
      dbdx = (2.0 * exp(-0.5 * x) * sinh(0.5 * x) - x) / (4.0 * sinh(0.5 * x)**2)
    else
      dbdx = -0.5 + x / 6.0
    end if
  end function

  elemental function phi1(x) result(phi)
    !! phi1(x) = (exp(x) - 1) / x = 1 / ber(x)
    real, intent(in) :: x
    real             :: phi

    phi = 1.0 / ber(x)
  end function

  elemental function dphi1dx(x) result(dphidx)
    !! derivative of phi1(x)
    real, intent(in) :: x
    real             :: dphidx

    dphidx = - dberdx(x) / (ber(x)**2)
  end function

  elemental function phi2(x) result(phi)
    !! phi2(x) = (exp(x) - 1 - x) / x**2
    real, intent(in) :: x
    real             :: phi

    ! local variables
    real,    parameter :: xr = 0.2
    integer, parameter :: m = 10
    integer            :: i
    real               :: fact

    if (abs(x) > xr) then
      phi = (exp(x) - 1 - x) / x**2
    else ! argument x is too small => use taylor series
      phi  = 0
      fact = 1
      do i = 0, m
        fact = fact * (i + 2)
        phi  = phi + x**i / fact
      end do
    end if
  end function

  elemental function dphi2dx(x) result(dphidx)
    !! derivative of phi2(x)
    real, intent(in) :: x
    real             :: dphidx

    ! local variables
    real,    parameter :: xr = 0.2
    integer, parameter :: m = 10
    integer            :: i
    real               :: fact

    if (abs(x) > xr) then
      dphidx = (exp(x) * (x - 2) + x + 2) / x**3
    else ! argument x is too small => use taylor series
      dphidx = 0
      fact   = 2
      do i = 1, m
        fact = fact * (i + 2)
        dphidx = dphidx + i * x**(i-1) / fact
      end do
    end if
  end function

  elemental function m_expm1(x) result(e)
    !! exp(x) - 1; accurate even for x close to 0
    real, intent(in) :: x
    real             :: e

    if (ieee_class(x) == IEEE_POSITIVE_INF) then
      e = x
    else
      e = exp(x)

      if (e == 1.0) then
        e = x
      else if (e - 1.0 == - 1.0) then
        e = -1
      else
        e = (e - 1.0) * x / log(e)
      end if
    end if
  end function

  elemental function m_log1p(x) result(l)
    !! log(1 + x); accurate even for x close to 0
    real, intent(in) :: x
    real             :: l

    real :: u, d

    if (ieee_class(x) == IEEE_POSITIVE_INF) then
      l = x
    else
      u = 1.0 + x
      d = u - 1.0

      if (d == 0) then
        l = x
      else
        l = log(u) * x / d
      end if
    end if
  end function

  pure function linspace(x0, x1, nx) result(x)
    !! create array of linear spaced values

    real,    intent(in) :: x0
      !! start value of x
    real,    intent(in) :: x1
      !! end value of x
    integer, intent(in) :: nx
      !! number of values
    real                :: x(nx)
      !! return array x

    integer :: i
    real    :: dx

    ! spacing between values
    dx = (x1 - x0) / (nx - 1)

    x(1) = x0
    do i = 2, nx - 1
      x(i) = x0 + (i - 1) * dx
    end do
    x(nx) = x1
  end function

  function logspace(x0, x1, nx) result(x)
    !! create array of logarithmic spaced values

    real,    intent(in) :: x0
      !! start value of x
    real,    intent(in) :: x1
      !! end value of x
    integer, intent(in) :: nx
      !! number of values
    real                :: x(nx)
      !! return array x

    integer :: i
    real    :: e0, e1, de

    ASSERT(x0 > 0.0)
    ASSERT(x1 > 0.0)

    e0 = log(x0)
    e1 = log(x1)
    de = (e1 - e0) / (nx - 1)

    x(1) = x0
    do i = 2, nx - 1
      x(i) = exp(e0 + (i - 1) * de)
    end do
    x(nx) = x1
  end function

  pure function eye_int(n) result(e)
    !! integer identity matrix

    integer, intent(in) :: n
      !! dimension
    integer             :: e(n,n)
      !! return integer identity matrix

    integer :: i

    e = 0
    do i = 1, n
      e(i,i) = 1
    end do
  end function

  pure function eye_real(n) result(e)
    !! real identity matrix

    integer, intent(in) :: n
      !! dimension
    real                :: e(n,n)
      !! return real identity matrix

    integer :: i

    e = 0.0
    do i = 1, n
      e(i,i) = 1.0
    end do
  end function

  pure function norm_inf(arr)
    !! Calculates the infinity norm $$ || \ \cdot \  ||_\infty $$ of the given array.

    real, intent(in) :: arr(:)
    real             :: norm_inf

    norm_inf = maxval(abs(arr))
  end function

  pure function check_lin_dep(M, rtol) result(l)
    !! Check if row vectors of M are linearly dependent
    real,           intent(in) :: M(:,:)
      !! input matrix
    real, optional, intent(in) :: rtol
      !! optional: relative tolerance (default: 1e-14)
    logical                    :: l
      !! return true if rows of M are lin. dep.; otherwise false

    integer :: nrows, ncols, ipiv(size(M,1)), row, row2, col, col2, pv, tmp
    real    :: rtol_, T(size(M,1),size(M,2)), nrm

    ! system size
    nrows = size(M,1)
    ncols = size(M,2)

    ! tolerance
    rtol_ = 1e-14
    if (present(rtol)) rtol_ = rtol

    ! init pivot indices
    do row = 1, nrows
      ipiv(row) = row
    end do

    ! work with copy of M
    T = M

    ! gauss elimination
    col = 0
    do row = 1, nrows
      col = col + 1

      ! find pivot
      do while (.true.)
        pv = row
        do row2 = row+1, nrows
          if (abs(T(ipiv(row2),col)) > abs(T(ipiv(pv),col))) pv = row2
        end do
        if (T(ipiv(pv),col) /= 0.0) exit
        col = col + 1
      end do

      ! virtually exchange current row with pivot
      if (pv /= row) then
        tmp       = ipiv(row)
        ipiv(row) = ipiv(pv)
        ipiv(pv ) = tmp
      end if

      ! gauss elimination for one column
      do row2 = row+1, nrows
        ! get row norm for zero check
        nrm = norm_inf(T(ipiv(row2),col:ncols))

        ! update row
        do col2 = col+1, ncols
          T(ipiv(row2),col2) = T(ipiv(row2),col2) - T(ipiv(row),col2) * T(ipiv(row2),col) / T(ipiv(row),col)
        end do

        ! check if row is zero => linearly dependent
        if (maxval(abs(T(ipiv(row2),col+1:ncols))) < rtol_ * nrm) then
          l = .true.
          return
        end if
      end do
    end do

    l = .false.
  end function

end module
